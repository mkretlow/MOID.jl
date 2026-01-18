# MOID.jl 

![status](https://img.shields.io/badge/status-deprecated-red)
[![Build Status](https://travis-ci.org/mkretlow/MOID.jl.svg?branch=master)](https://travis-ci.org/mkretlow/MOID.jl)

This package has been **renamed and superseded by AstroMOID.jl**.

Please use:
https://github.com/mkretlow/AstroMOID.jl

This repository is kept for archival reasons only.
No further development or bug fixes will happen here.

___

Compute the MOID - Minimum Orbit Intersection Distance for two given confocal, elliptical orbits.
It uses the idea of rotating meridional plane and calculates the MOIDs numerically.

The method is described in the paper:

T.Wiśniowski and H.Rickman, "A Fast, Geometric Method for Calculating Accurate Minimum Orbit Intersection Distances (MOIDs)", Acta Astronomica, Vol. 63 (2013) pp. 293–307. The paper and the original Fortran source code are provided in the `references/` subdir.


## Install

```julia
using Pkg
Pkg.add("MOID")
```

Or in the Julia REPL package mode (press `]`):

```julia
pkg> add MOID
```

## Quickstart

Calculate the MOID between asteroids (1) Ceres and (30) Urania. Argument values are:

- semi-major axis (a in au)
- eccentricity (e)
- argument of perhielion (ω)
- longitude of ascending node (Ω)
- inclination (i) (all angles in rad)

for both objects. Result of the `wisric_moid` function is the MOID in au.

```julia
using MOID

rad = pi/180
ceres = [2.7691652, 0.0760091, 73.59764*rad, 80.30553*rad, 10.59407*rad]
urania = [2.3655722, 0.127581, 87.42605*rad, 307.46872*rad, 2.09575*rad]

wisric_moid(ceres...,urania...)

# Output: 0.24521440655831886
```

---
