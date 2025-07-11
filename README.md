# Propagator created: Wed 10. July 2025, 17:05 PM

A time-propagation using C(t+dt) = \exp(-i dt H(t)) C(t)
## Prerequisites

- CMake ≥ 3.10
- A C/C++ compiler
- Make

## Build

```bash
# 1. Create and enter build directory
mkdir -p build
cd build

# 2. Configure
cmake ..

# 3. Compile (use all cores)
make -j$(nproc)
