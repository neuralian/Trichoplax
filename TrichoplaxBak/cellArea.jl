
function cellArea(x::Array{Float64,1}, y::Array{Float64,1}, n::Int64)
# area of n-gon with vertices x,y
# vertices must be ordered - clockwise or anticlockwise

    A = 0.0
    for i in 1:n
        j = i % 12 + 1
        A = A + x[i]*y[j] - y[i]*x[j]
    end

    return (abs(A/2.0))

end
