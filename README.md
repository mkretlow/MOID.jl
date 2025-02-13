# MOID.jl [![Build Status](https://travis-ci.org/mkretlow/MOID.jl.svg?branch=master)](https://travis-ci.org/mkretlow/MOID.jl)
Compute the MOID - Minimum Orbit Intersection Distance for two given confocal, elliptical orbits.
It uses the idea of rotating meridional plane and calculates the MOIDs numerically.

The method is described in the paper:

T.Wiśniowski and H.Rickman, "A Fast, Geometric Method for Calculating Accurate Minimum Orbit Intersection Distances (MOIDs)", Acta Astronomica, Vol. 63 (2013) pp. 293–307. The paper and the original Fortran source code are provided here in the docs subdir.


## Install

```julia
pkg> add "https://github.com/mkretlow/MOID.jl.git"
```

## Quickstart
```julia
julia> using MOID

# Calculate the MOID between asteroids (1) Ceres and (30) Urania. Argument values are: 
# semi-major axis (au), eccentricity, argument of perhielion (ω), longitude of ascending node (Ω),
# inclination (i) (all angles in rad). Result is MOID in au.

julia> rad = pi/180
julia> ceres = [2.7691652, 0.0760091, 73.59764*rad, 80.30553*rad, 10.59407*rad]
julia> urania = [2.3655722, 0.127581, 87.42605*rad, 307.46872*rad, 2.09575*rad]
julia> wisric_moid(ceres...,urania...)
0.24521440655831886
```

---
