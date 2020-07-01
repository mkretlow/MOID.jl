# MOID.jl [![Build Status](https://travis-ci.com/mkretlow/MOID.jl.svg?branch=master)](https://travis-ci.com/mkretlow/MOID.jl)
Compute the MOID - Minimum Orbit Intersection Distance for two given confocal, elliptical orbits.
It uses the idea of rotating meridional plane and calculates the MOIDs numerically.

This Julia module is a wrapper to the Fortran program by the original authors Hans Rickman  and Tomasz Wiśniowski, both
Space Research Centre, Polish Academy of Sciences.

The method is described in their paper:

T.Wiśniowski and H.Rickman, “A Fast, Geometric Method for Calculating Accurate Minimum Orbit Intersection Distances (MOIDs)” 2013 Acta Astronomica

The paper and the source code are provided here in the deps and docs subdirs.


## Install
MOID.jl can then be installed through Julia's package manager. On Linux and macOS you need to have either Gfortran or the Intel Fortran Compiler installed to be able to build the binary dependencies.

```julia
pkg> add "https://github.com/mkretlow/MOID.jl.git"
```

## Quickstart
```julia
julia> using MOID

# Calc MOID between asteroids (1) Ceres and (30) Urania. Argument values are: 
# semi-major axis (AU), eccentricity, argument of perhielion (ω), longitude of ascending node (Ω),
# inclination (i) (all in degrees). Result is MOID in AU.

julia> ceres = [2.7691652, 0.0760091,  73.59764,  80.30553,10.59407]
julia> urania = [2.3655722, 0.127581 ,  87.42605, 307.46872, 2.09575]
julia> moidF(ceres...,urania...)
0.24521440655831864
```

---
