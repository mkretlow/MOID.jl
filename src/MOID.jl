# ==================================================================================================================================
# Module https://github.com/mkretlow/MOID.jl : Mike Kretlow [astrodynamics.de], Start 2018, MIT "Expat" License
# ==================================================================================================================================

module MOID

using Libdl
using Printf

export calc_moidF


# Load shared library

const lib = find_library(["libmoid"],["deps", joinpath(dirname(@__FILE__), "..", "deps")])


if isempty(lib)
	error("Could not find shared library")
end

@info("Using shared library $lib")



"""
   `moid = calc_moidF(sma1, ecc1, peri1, node1, incl1, sma2, ecc2, peri2, node2, incl2)`

Calculating Accurate Minimum Orbit Intersection Distances (MOID) by calling original Fortran routine

## Args

   sma1, ecc1, peri1, node1, incl1, sma2, ecc2, peri2, node2, incl2 [angles in deg !]

## Return
   MOID in units like sma (usually AU)

## Reference

T.Wisniowski and H.Rickman "A Fast, Geometric Method for Calculating Accurate Minimum Orbit Intersection Distances (MOIDs)"
published in 2013 in Acta Astronomica.

See also: http://moid.cbk.waw.pl/orbity/default/index

"""
function calc_moidF(sma1::Float64, ecc1::Float64, peri1::Float64, node1::Float64, incl1::Float64,
                     sma2::Float64, ecc2::Float64, peri2::Float64, node2::Float64, incl2::Float64)

   ccall((:moid_,lib),Float64,(Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},
                               Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64}),
                               sma1, ecc1, peri1, node1, incl1, sma2, ecc2, peri2, node2, incl2)
end

end  # module


#==================================================================================================================================#