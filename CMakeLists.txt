cmake_minimum_required(VERSION 3.13)
project(propagator LANGUAGES CXX)

# ── User-tweakable options ──────────────────────────────────────
set(CMAKE_CXX_STANDARD            17)
set(CMAKE_CXX_STANDARD_REQUIRED   ON)
set(CMAKE_CXX_EXTENSIONS          OFF)

# Change to "Release" if you always want optimisation
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# ── Source files ────────────────────────────────────────────────
add_executable(propagate
    propagate.cpp      # main program (includes propagate.hpp)
)

# Optimise and warn reasonably
target_compile_options(propagate
    PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-O2 -Wall -Wextra -Wpedantic>
        $<$<CXX_COMPILER_ID:MSVC>:/O2 /W4>
)

# ── Installation (optional) ─────────────────────────────────────
#   make install   → installs the executable to CMAKE_INSTALL_PREFIX/bin
install(TARGETS propagate
        RUNTIME DESTINATION bin)

# ── Packaging stub (optional) ───────────────────────────────────
#   cpack           → produces a .tar.gz or .zip in the build dir
include(CPack)
