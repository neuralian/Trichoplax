# Trichoplax development model
# MGP Dec 2019

using Makie
using LinearAlgebra
using Statistics

struct Trichoplax
    x0::Float64       # un-stressed edge length
    v0::Float64       # un-stressed cell volume (area in 2D)
    vertex::Array    # nVertex*2 cell vertices
    cell::Array      # nCells*6  indices of cell vertices
    edge::Array      # nLink*2   indices of 2 cell vertices
end

diameter = 100

function make_trichoplax(Diameter)
    # construct a Trichoplax of specified diameter

    plaxDiam = Float64(Diameter)  # Trichoplax diameter microns
    cellDiam = 5.0   # cell diameter microns

    # row height factor
    # (row height is h*cellDiam, column width is cellDiam)
    h = sqrt(3.0)/2.0

    # vertex distance tolerance
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
    cellIndex= fill(-1, MaxNCells, 6)

    # links
    # each row defines a link between cell vertices
    # by specifying a pair of row indices in the vertex array
    vertexLink = fill(-1, 6*MaxNCells, 2)

    # cell nuclei
    # used to construct cell membranes/cells, not members of Trichoplax
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
                    cellIndex[cell, i] = j
                    vertexExists = true
                end
            end
            if !vertexExists
                nCellVertex = nCellVertex + 1
                cellVertex[nCellVertex, :] = newVertex
                cellIndex[cell, i] = nCellVertex
            end
            # create links
            if i>1
                linkExists = false
                for j in 1:nLink
                    if (transpose([cellIndex[cell,i-1] cellIndex[cell,i]]') ==
                                                          vertexLink[j,:]') |
                       (transpose([cellIndex[cell,i] cellIndex[cell,i-1]]') ==
                                                          vertexLink[j,:]')
                        linkExists = true
                    end
                end
                if !linkExists
                    nLink = nLink + 1
                    vertexLink[nLink,:] = [cellIndex[cell,i-1] cellIndex[cell,i]]
                end
            end
            if nLink>0
                lines!(cellVertex[vertexLink[nLink, :],1],
                    cellVertex[vertexLink[nLink, :],2])
                # display(petridish)
                # sleep(1/25)
            end
       end
       # add 6th edge (close the loop from 6th to 1st vertex)
       linkExists = false
       for j in 1:nLink
           if (transpose([cellIndex[cell,1] cellIndex[cell,6]]') ==
                                                       vertexLink[j,:]') |
              (transpose([cellIndex[cell,6] cellIndex[cell,1]]') ==
                                                        vertexLink[j,:]')
               linkExists = true
           end
       end
       if !linkExists
           nLink = nLink + 1
           vertexLink[nLink,:] = [cellIndex[cell,1] cellIndex[cell,6]]
       end
    end

    trichoplax = Trichoplax(h, sqrt(3.)/2.0*cellDiam^2,
                            cellVertex[1:nCellVertex,:],
                            cellIndex[1:nCells,:],
                            vertexLink[1:nLink,:])

end


function draw_trichoplax(trichoplax, scene)

    nVertex = size(trichoplax.vertex, 1)
    nCell = size(trichoplax.cell, 1)
    vertex = fill(0.0, nVertex+nCell, 2) # will add centre vertex to each cell
    for i in 1:nVertex
        vertex[i,:] = trichoplax.vertex[i,:]
    end
    facet = fill(0,6*nCell,3) # 6 triangle facets per cell
    nFacet = 0
    for cell in 1:nCell
        vertex[nVertex+cell, :] = mean(vertex[trichoplax.cell[cell,:],:], dims=1)
        for i in 1:6
            nFacet = nFacet + 1
            i0 = trichoplax.cell[cell,i]
            i1 = trichoplax.cell[cell,i%6+1]
            facet[nFacet, :] = [nVertex+cell i0 i1]
        end
    end



            poly!(vertex, facet, color =  rand(nVertex+nCell) )



    for link in 1:size(trichoplax.edge,1)
        lines!(scene, trichoplax.vertex[trichoplax.edge[link, :],1],
               trichoplax.vertex[trichoplax.edge[link, :],2])
    end



    display(scene)
end



trichoplax = make_trichoplax(diameter);


# Set scene
sceneWidth = Int64(round(diameter*1.5))
if isodd(sceneWidth)
    sceneWidth = sceneWidth + 1
end

petridish = Scene(limits = FRect(-sceneWidth/2, -sceneWidth/2,
                                sceneWidth,sceneWidth), scale_plot = false)

draw_trichoplax(trichoplax, petridish)


# # enclosing circle
# nPt = 128.
# lines!(diameter/2.0*cos.(2π*(0:nPt)./nPt),
#         diameter/2.0*sin.(2π*(0:nPt)./nPt),
#         color = :lightblue)




# println(MaxNCells)

function cellVolume(cell, trichoplax)
# volume (area in 2D) of cell

    V = 0.0
    for i in 1:6
        j = i % 6 + 1
        V = V + trichoplax.vertex[trichoplax.cell[cell,i], 1] *
                trichoplax.vertex[trichoplax.cell[cell,j], 2] -
                trichoplax.vertex[trichoplax.cell[cell,i], 2] *
                trichoplax.vertex[trichoplax.cell[cell,j], 1]
    #    A = A + x[i]*y[j] - y[i]*x[j]
    end

    return (abs(V/2.0))

end

function cellVolume2(cell, trichoplax)
# volume (area in 2D) of cell

    V = 0.0
    for i in 1:5
        V = V + trichoplax.vertex[trichoplax.cell[cell,i], 1] *
                trichoplax.vertex[trichoplax.cell[cell,i+1], 2] -
                trichoplax.vertex[trichoplax.cell[cell,i], 2] *
                trichoplax.vertex[trichoplax.cell[cell,i+1], 1]
    end

    V = V + trichoplax.vertex[trichoplax.cell[cell,6], 1] *
            trichoplax.vertex[trichoplax.cell[cell,1], 2] -
            trichoplax.vertex[trichoplax.cell[cell,6], 2] *
            trichoplax.vertex[trichoplax.cell[cell,1], 1]

    return (abs(V/2.0))

end


function stress(v, trichoplax)
    # mechanical stress on trichoplax with vertices v

    # turgor pressure
    E = 0.0
    for i in 1:size(v,1)
        E = E + (cellVolume(:, trichoplax)-trichoplax.v0)^2
    end
end
