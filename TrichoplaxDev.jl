# Trichoplax development model
# MGP Dec 2019

using Makie
using LinearAlgebra

function trichoplax(Diameter)

plaxDiam = Float64(Diameter)  # Trichoplax diameter microns
cellDiam = 5.0   # cell diameter microns

# row height factor
# (row height is h*cellDiam, column width is cellDiam)
h = sqrt(3.0)/2.0

# vertex tolerance (vertices closer than this are combined into one)
vertexTol = cellDiam/100.

nRows     = Int64(ceil(plaxDiam/(cellDiam*h)/2.0))
nCols     = Int64(ceil(plaxDiam/cellDiam/2.0))
MaxNCells = Int64(ceil(π*(nRows+1)*(nCols+1)))

# cell centre x-y coordinate array
cellNucleus = NaN*ones(MaxNCells,2)

# cell vertex array
# each row is a vertex (x,y) of a hexagonal cell boundary
# TODO: Determine empirically how big (or small) this array needs to be
cellVertex = NaN*ones(6*MaxNCells, 2)

# cell boundary indices
# each hexagonal cell defined by 6 vertex indices listed anticlockwise
cellFace = fill(-1, MaxNCells, 6)

# links
# each row defines a link between cell vertices
# by specifying a pair of row indices in the vertex array
vertexLink = fill(-1, 6*MaxNCells, 2)

# Set scene
sceneWidth = Int64(round(plaxDiam*1.5))
if isodd(sceneWidth)
    sceneWidth = sceneWidth + 1
end
petridish = Scene(limits = FRect(-sceneWidth/2, -sceneWidth/2,
                                sceneWidth,sceneWidth), scale_plot = false)

# draw cell nuclei
nCells = 0
for row in 0:nRows
    for col in 0:nCols
        if iseven(row)
            x = col*cellDiam
            y = h*row*cellDiam
            if (x^2 + y^2) <= (plaxDiam/2.)^2
                nCells = nCells + 1
                cellNucleus[nCells,:] = [x y]
                if y>cellDiam/4.
                    nCells = nCells + 1
                    cellNucleus[nCells,:] = [x -y]
                end
                if x > cellDiam/4.
                    nCells = nCells + 1
                    cellNucleus[nCells,:] = [-x y]
                end
                    if (x > cellDiam/4.) &  (y > cellDiam/4.)
                    nCells = nCells + 1
                    cellNucleus[nCells,:] = [-x -y]
                end
            end
        else
            x = (col + 0.5)*cellDiam
            y = h*row*cellDiam
            if (x^2 + y^2) <= (plaxDiam/2.)^2
                nCells = nCells + 1
                cellNucleus[nCells,:] = [x y]
                nCells = nCells + 1
                cellNucleus[nCells,:] = [-x -y]
                nCells = nCells + 1
                cellNucleus[nCells,:] = [x -y]
                nCells = nCells + 1
                cellNucleus[nCells,:] = [-x y]
            end
        end
    end
end

# draw cell centres
scatter!(petridish, eachcol(cellNucleus[1:nCells,:])...,
                 markersize = cellDiam/8., color = :grey)

# construct cells and links
nCellVertex = 0
nLink = 0
for cell in 1:nCells
    for i in 1:6
        newVertex = cellNucleus[cell,:]' .+
                     0.5/h*cellDiam.*[cos(2π*(i-.5)/6)  sin(2π*(i-.5)/6)]
        # check for existing vertex
        vertexExists = false
        for j in 1:nCellVertex
            if norm(cellVertex[j,:]' - newVertex) < vertexTol
                cellFace[cell, i] = j
                vertexExists = true
            end
        end
        if !vertexExists
            nCellVertex = nCellVertex + 1
            cellVertex[nCellVertex, :] = newVertex
            cellFace[cell, i] = nCellVertex
        end
        # create links
        if i>1
            linkExists = false
            for j in 1:nLink
                if [cellFace[cell,i-1] cellFace[cell,i]] == vertexLink[j,:] or
                [cellFace[cell,i] cellFace[cell,i-1]] == vertexLink[j,:]
                    linkExists = true
                end
            end
            if !linkExists
                nLink = nLink + 1
                vertexLink[nLink,:] = [cellFace[cell,i-1] cellFace[cell,i]]
            end
        end
   end
end

for i in 1:nCells
    lines!(cellVertex[vcat(cellFace[i,:],cellFace[i,1]),1],
            cellVertex[vcat(cellFace[i,:],cellFace[i,1]),2])
end



# enclosing circle
nPt = 128.
lines!(plaxDiam/2.0*cos.(2π*(0:nPt)./nPt),
        plaxDiam/2.0*sin.(2π*(0:nPt)./nPt),
        color = :lightblue)


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
# println(nCells)
# println(MaxNCells)

end
