# Draw hexagon as N-layer hexagonal grid
# (C) Mike Paulin University of Otago 2019

using Printf
using Makie

N = 64

NN = Int64(round(N^2 - N*(N-1)/4))
nuc = fill(0.0,NN, 2)

h = sqrt(3.)/2.
x = 0.
y = 0.
cellDiam = 5.0
s = Scene(scale_plot = false, markersize = cellDiam/4)
nCells = 0
for i in 1:N
    global y = cellDiam*(i-1)*(1+h)/2.
    for j in 1:Int64(ceil((2N-i)/2))
        global x = cellDiam*(j - 0.5 - 0.5*isodd(i))
    #@printf("%s%.1f%s%.3f%s", "(", x, ",", y, ")")
    #Makie.scatter!([x],[y])
    global nCells = nCells+1
    global nuc[nCells, :] = [x,y]
    end
    #println("")
end
Makie.scatter!(nuc[:,1], nuc[:,2])
display(s)
println((nCells, NN))
