# Trichoplax development model
# MGP Dec 2019

using Makie


MaxNumCells = 100
cellRadius = 0.5

# declare array for x-y coordinates of cell centres
nucleus = NaN*ones(MaxNumCells,2)

# declare array for neighbour links
# each cell links to up to 6 neighbours
neighbour = fill(-1,MaxNumCells, 6)

# declare array for x-y coordinates of cell vertices
# (cell-cell binding locations = desmosomes)
# these are centroids of triangular faces formed by
#    neighbour links, plus some extras to enclose
#    nuclei around the perimeter.  The number of extras
#    required will be determined empirically. In the
#    meantime we use an upper bound
nEdgeVertices = MaxNumCells
desmosome = NaN*ones(MaxNumCells+nEdgeVertices, 2)

# declare array for cell membranes
# (Voronoi tesselation)
# each cell has six edges, defined by 6 points anticlockwise
cell = NaN*ones(MaxNumCells, 6)


# initialize with 3 cells (3 nuclei = 1 Delaunay triangle)
nCells = 3
for i in 1:3
    a = π*(1.0/2.0 + 2.0*(i-1)/3.0)
    nucleus[i,:] = cellRadius.*[cos(a), sin(a)]
    neighbour[i,1:2] = [i i+1].%3 .+1
end

# desmosome at origin is the common vertex of the three initial cells
desmosome[1,:] = [0.0 0.0]




# draw initial Trichoplax
scene = Scene(limits = FRect(-5.0, -5.0, 10.0, 10.0),
              scale_plot = false)

# draw cell nuclei
scatter!(scene, eachcol(nucleus[1:nCells,:])...,
                markersize = 0.2, color = :grey)
# draw neighbour links
# nb this draws two links on top of each other for each pair of neighbours
#    but it's not worth the hassle to draw only one
for i in 1:nCells
    for j in 1:6
        if neighbour[i,j]>0    # ith cell has a jth neighbour
            lines!(scene,
              [nucleus[i,1], nucleus[neighbour[i,j], 1]],
              [nucleus[i,2], nucleus[neighbour[i,j], 2]],
              color = :grey)
        end
    end
end




display(scene)
