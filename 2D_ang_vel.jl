using CairoMakie
using LinearAlgebra
using Statistics
using SpecialFunctions

# the following code takes in the number of grid points, bdys and returns a matrix with the non-dimensional u
# at each grid point. The formula is taken from Charbonneau (2010) as mentioned in the report

function TwoD_diff_rot(mp, np, x_bdy, y_bdy, invar; Ωc = 0.8752, a2 = 0.1264, a4 = 0.1591, w = 0.1)

    # in spherical it depends on r and θ, so here it depends on z and x

    output = zeros(mp,np)
    x_ra = LinRange(x_bdy[1], x_bdy[2], np)
    y_ra = LinRange(y_bdy[1], y_bdy[2], mp)

    for yc in 1:mp
        for xc in 1:np
            
            output[yc,xc] = Ωc + (1 + erf(y_ra[yc]/w))*( (1-a2*(cos(x_ra[xc]))^2 - a4*(cos(x_ra[xc])^2) ) - Ωc )/(2) 

        end
    end

    return output

end
