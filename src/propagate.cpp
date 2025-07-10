//propagate.cpp
#include <complex>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <chrono>
#include "propagate.hpp"

using cxd   = std::complex<double>;
//–––––––––––––– helpers ––––––––––––––
double l2norm(const std::vector<cxd>& v){
    double s = 0.0;
    for(int i = 0; i < v.size(); i++) s+= std::norm(v[i]);
    return std::sqrt(s);
}

//–––––––– physical constants
constexpr double au_time_fs = 0.02418884326505;// 1 a.u. time  in fs
constexpr double I0_Wpcm2   = 3.50944506e16;

//–––––––– pulse shape
double laser_pulse(double t_au,double I_Wpcm2, double T_au, double omega){
    double F0 = std::sqrt(I_Wpcm2/I0_Wpcm2);
    double A0 = F0/omega;
    double x = t_au/T_au;
    return A0*std::exp(-x*x)*std::cos(omega*t_au);
}

//–––––––– H·C
std::vector<cxd> H_times_C(const std::vector<double>& E,
                const std::vector<std::vector<double> >& D,const std::vector<cxd>& C,
                double t,double I_Wpcm2,double T, double omega){
    const size_t N = E.size();
    std::vector<cxd> result(N, cxd{0.0, 0.0});
    double A = laser_pulse(t, I_Wpcm2, T, omega);

    for(size_t j = 0; j < N; ++j){
        result[j] += C[j] * E[j]; // diagonal
        for(size_t k = 0; k < N; ++k){
            if(D[j][k] == 0.0) continue;
            cxd Vjk = cxd{0.0, (E[j] - E[k]) * D[j][k] * A};
            result[j] += Vjk * C[k];
        }
    }
    return result;
}

//–––––––– single Taylor step
std::vector<cxd> propagate_step(const std::vector<double>& E,
                    const std::vector<std::vector<double> >& D,
                    const std::vector<cxd>&   C,
                    double                    t,
                    double                    dt,
                    double                    I_Wpcm2,
                    double                    T,
                    double                    omega,
                    double                    tol = 1e-12,
                    int                       max_order = 100)
{
    cxd minus_i_dt{0.0, -dt};
    std::vector<cxd> Cnext = C;

    std::vector<cxd> HtC = H_times_C(E, D, C, t, I_Wpcm2, T, omega);
    std::vector<cxd> term(HtC.size(), cxd{0.0, 0.0});
    for(int k = 0; k < HtC.size(); k++){
        term[k] = minus_i_dt*HtC[k];
    }

    int n = 1;
    while(l2norm(term) > tol * l2norm(Cnext) && n < max_order)
    {
        for(int i = 0; i < Cnext.size(); ++i) Cnext[i] += term[i];
        std::vector<cxd> HtC_t = H_times_C(E, D, term, t, I_Wpcm2, T, omega);
        for(int k = 0; k < HtC_t.size(); k++){
            term[k] = ( minus_i_dt/double(n + 1) )*HtC_t[k];
        }
        ++n;
    }
    for(int i = 0; i < Cnext.size(); ++i) Cnext[i] += term[i];
    //std::cout<<"Taylor expansion terminated at order n = "<<n<<std::endl;
    return Cnext;
}

//–––––––– full propagation loop
std::vector<cxd> run_propagation(const std::vector<double>& E,
                     const std::vector<std::vector<double> >& D,
                     double                    I_Wpcm2,
                     double                    omega,
                     double                    T_au,
                     double                    dt,
                     double                    t_start,
                     double                    t_final)
{
    const int N = E.size();
    std::vector<cxd> C(N, cxd{0.0, 0.0});
    C[0] = 1.0; // ground state; C[1]=1.0 is excited state

    std::cout << std::fixed << std::setprecision(10);
    int step = 0;
    double t = t_start;
    while(t < t_final){
        if(step % 100 == 0)
            std::cout<<"t = "<<t*au_time_fs<<" fs   |C|² = " <<l2norm(C)*l2norm(C)<<std::endl;

        std::vector<cxd> Cnext = propagate_step(E, D, C, t, dt, I_Wpcm2, T_au, omega);

        double norm = l2norm(Cnext);
        if(std::abs(norm - 1.0) > 1e-6) std::cerr<<"Warning: |C| drift = "<<norm<<" at t = "<<t<<std::endl;

        C.swap(Cnext);//C = std::move(Cnext);
        t += dt;
        ++step;
    }
    return C;
}

//–––––––––––––––––––– main
int main(){

    std::cout<<"1. initialization...\n";

    auto t_start = std::chrono::high_resolution_clock::now();
    //locations & filenames
    const std::string eigen_file  = "../helium/eigenvalues_2000_2_5.dat";
    const std::string dipole_file = "../helium/Integraldipole_2000_2_5.dat";

    //read data
    int l_max = 0;
    size_t len_E = 0;
    auto E = read_energies(eigen_file, l_max, len_E);
    auto D = read_dipoles(dipole_file, E.size());

    //laser / propagation parameters
    const double I_Wpcm2 = 1.0e12;                     //intensity
    const double omega   = 1.5;                       //freq. (a.u.)
    const double T_au    = 2.0/au_time_fs;           //3 fs
    const double dt      = 0.1;                     //time step
    const double t0      = -3.0 * T_au;            //start time
    const double tf      = 3.0 * T_au;            //stop time

    auto t_end_ini = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_ini = t_end_ini - t_start;
    std::cout<<"              took "<<elapsed_ini.count()/60.0<<" mins.\n";

    auto t_start_exe = std::chrono::high_resolution_clock::now();
    std::cout<<"2. execution...\n";
    std::vector<cxd> C_final = run_propagation(E, D,I_Wpcm2, omega, T_au,dt, t0, tf);

    auto t_end_exe = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_exe = t_end_exe - t_start_exe;
    std::cout<<"Execution took "<<elapsed_exe.count()/60.0<<" mins.\n";

    std::cout<<"3. post-processing...\n";
    //Write the spectrum
    spectrum::write_spectrum(E,len_E,l_max,C_final,/*e_up*/E.back(),/*grid spacing dE*/0.0025); 

    auto t_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = t_end - t_start;
    std::cout<<"4. done! running took "<<elapsed.count()/60.0<<" mins.\n";


    return 0;
}
