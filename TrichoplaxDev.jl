# Trichoplax development model
# MGP Dec 2019

using Makie

plaxDiam = 100.0  # Trichoplax diameter in cell diameters
cellDiam = 5.0   # cell diameter microns

# row height factor
# (row height is h*cellDiam, column width is cellDiam)
h = sqrt(3.0)/2.0

nRows     = Int64(ceil(plaxDiam/(cellDiam*h)))
nCols     = Int64(ceil(plaxDiam/cellDiam))
MaxNCells = Int64(ceil(π*(nRows+1)*(nCols+1)))

# array for x-y coordinates of cell centres
cellNucleus = NaN*ones(MaxNCells,2)

# cell vertex array
# each row is a vertex (x,y) of a hexagonal cell boundary
# TODO: Determine empirically how big (or small) this array needs to be
cellVertex = NaN*ones(25*MaxNCells, 2)

# array for cell boundaries
# each hexagonal cell defined by 6 vertex indices listed anticlockwise
cellBoundary = fill(1, MaxNCells, 6)

# Set scene
petridish = Scene(limits = FRect(-2*plaxDiam, -2*plaxDiam, 4*plaxDiam,4*plaxDiam),
              scale_plot = false)

# draw cell nuclei
nCells = 1
cellVertex[1,:] = [0.0 0.0]
for row in 0:nRows
    for col in 0:nCols
        if iseven(row)
            x = col*cellDiam
            y = h*row*cellDiam
            if (x^2 + y^2) < (plaxDiam + cellDiam)^2
                global nCells = nCells + 1
                cellVertex[nCells,:] = [x y]
                if row>0
                    nCells = nCells + 1
                    cellVertex[nCells,:] = [x -y]
                end
                if col>0
                    nCells = nCells + 1
                    cellVertex[nCells,:] = [-x y]
                end
                if (row>0) & (col>0)
                    nCells = nCells + 1
                    cellVertex[nCells,:] = [-x -y]
                end
            end
        else
            x = (col - 0.5)*cellDiam
            y = h*row*cellDiam
            if (x^2 + y^2) < (plaxDiam + h*cellDiam)^2
                nCells = nCells + 1
                cellVertex[nCells,:] = [x y]
                if row>0
                    nCells = nCells + 1
                    cellVertex[nCells,:] = [x -y]
                end
                if col>0
                    nCells = nCells + 1
                    cellVertex[nCells,:] = [-x y]
                end
                if (row>0) & (col>0)
                    nCells = nCells + 1
                    cellVertex[nCells,:] = [-x -y]
                end
            end
        end
    end
end


scatter!(petridish, eachcol(cellVertex[1:nCells,:])...,
                 markersize = 1, color = :grey)


# # draw neighbour links
# # nb this draws both links on top of each other for each pair of neighbours
# #    but it's not worth the hassle to draw only one
# for i in 1:nCells
#     for j in 1:6
#         if neighbour[i,j]>0    # ith cell has a jth neighbour
#             lines!(scene,
#               [nucleus[i,1], nucleus[neighbour[i,j], 1]],
#               [nucleus[i,2], nucleus[neighbour[i,j], 2]],
#               color = :grey)
#         end
#     end
# end
#
# # vertex at origin is the common vertex of the three initial cells
# desmosome[1,:] = [0.0 0.0]
#
# # Cell membrane vertex list for each cell (nucleus)
# # The initial three cells initially share the first vertex
# #cell[1:3, :] = [1 0 0 0 0 0]
#
# # construct cell vertices
# for i in 1:nCells
#
#
#
# end


display(petridish)
println(nCells)
println(MaxNCells)
