//propagate.hpp
#ifndef PROPAGATE_HPP
#define PROPAGATE_HPP

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <array>

//––––––––––––––––––––––––– FILE READERS ––––––––––––––––––––––––––
// 3-column eigenvalue file  l  n  E
std::vector<double>read_energies(const std::string& file, int& l_max_out, std::size_t& len_E_out)
{
    std::ifstream in(file);
    if(!in) throw std::runtime_error("Cannot open " + file);

    std::vector<std::vector<double>> by_l;
    double l_in, n_in, E_in;
    while(in >> l_in >> n_in >> E_in){
        int l = static_cast<int>(l_in);
        if(l >= static_cast<int>(by_l.size()))
            by_l.resize(l + 1);
        by_l[l].push_back(E_in);
    }

    if(by_l.empty()) throw std::runtime_error("No data in " + file);

    len_E_out = by_l[0].size();
    l_max_out = static_cast<int>( by_l.size() );

    for(const auto& v : by_l) if(v.size() != len_E_out) throw std::runtime_error("Different n-count per l in " + file);

    std::vector<double> E;
    E.reserve(static_cast<std::size_t>(l_max_out) * len_E_out);
    for(int l = 0; l < l_max_out; ++l)
        for(double e : by_l[l])
            E.push_back(e);
    return E;                     // flat (l-major)
}


// Flat N×N dipole matrix
std::vector<std::vector<double> > read_dipoles(const std::string& file, std::size_t N)
{
    std::vector<std::vector<double> > D(N, std::vector<double>(N));
    std::ifstream din(file);
    if(!din)throw std::runtime_error("Cannot open " + file);

    for(std::size_t i = 0; i < N; ++i){
        for(std::size_t j = 0; j < N; ++j){
            if( !(din >> D[i][j]) ){
                throw std::runtime_error("Not enough entries in " + file);
            }
        }
    }
    return D;
}


//-----------------------------------------------------------------
//  Photo-electron spectrum  (Gaussian broadening)
//-----------------------------------------------------------------
namespace spectrum{

    constexpr double pi = 3.14159265358979323846;

    /// Gaussian kernel  K(E; dE)  normalised to ∫K dE = 1
    inline double K_Gaussian(double E, double dE){
        const double norm = 1.0/(std::sqrt(pi) * std::abs(dE));
        return norm*std::exp(-(E * E)/(dE * dE));
    }

    //build  occ[l][n]  from the flattened C
    std::vector<std::vector<double> >
    populations_from_C(const std::vector<std::complex<double> >& C,int l_max,std::size_t len_E){
        std::vector<std::vector<double> > occ(l_max,std::vector<double>(len_E, 0.0));
        for(int l = 0; l < l_max; ++l)
            for(std::size_t n = 0; n < len_E; ++n)
                occ[l][n] = std::norm(C[l * len_E + n]);
        return occ;
    }
    /// Write photoelectron spectrum
    /// • E        – flattened energies  E[l*len_E + n]  (a.u.)
    /// • len_E    – number of box states per l
    /// • l_max    – number of l-channels actually present
    /// • occ      – |C|² populations  occ[l][n]  (same shape as E)
    /// • e_up     – output energy ceiling (a.u.)
    /// • step     – energy grid spacing (a.u.)
    ///
    /// The file has seven columns:
    ///   Es,  P_tot,  P_l=0,  P_l=1,  … , P_l=4
    /// and is named  spectrum.dat
    ///
    void write_spectrum(const std::vector<double>& E,
                        std::size_t len_E, int l_max,
                        const std::vector<std::complex<double> >& C_final,
                        double e_up, double step= 0.0025){
        //1.Convert amplitudes → populations
        auto occ = populations_from_C(C_final, l_max, len_E);
        std::ofstream pout("spectrum.dat");
        pout<<std::fixed << std::setprecision(15);

        for(double Es = step; Es < e_up; Es += step)
        {
            double P_tot     = 0.0;
            std::array<double,5> P_l{};           // l = 0…4 (auto-zeroed)

            for(int l = 0; l < l_max; ++l)
            {
                for(std::size_t n = 0; n + 1 < len_E; ++n){
                    const double En = E[l * len_E + n];
                    if (En <= 0.0) continue;       // bound states → ignore

                    const double dE = E[l * len_E + n + 1] - En;   // local spacing
                    const double w  = occ[l][n]*K_Gaussian(Es - En, dE);

                    P_tot+= w;
                    if(l < 5) P_l[l] += w;       // up to l = 4 for extra columns
                }
            }

            pout<<Es<<'\t'<<P_tot;
            for(double val : P_l)pout<<'\t'<<val;
            pout<< '\n';
        }
    }
} // namespace spectrum


#endif // PROPAGATE_HPP
