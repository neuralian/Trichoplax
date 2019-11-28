using ForwardDiff

function areaXY(x)
    # area of polygon whose vertices are specified in a column vector
    # containing all the x-coords followed by all the y-coords.
    # This function is intended to make it easy to do calculus (numerical
    # optimization) of functions of the area of a polygon (e.g. surface energy).
    # NB for speed it does no argument checking: caveat usor
    # MGP Nov 2019

    A = 0.0
    n = Int64(length(x)/2)
    for i in 1:(n-1)
         A = A + x[i]*x[n+i+1] - x[i+1]*x[n+i]
    end
    A = A + x[n]*x[n+1] - x[2*n]*x[1]

    return A/2.0

end
