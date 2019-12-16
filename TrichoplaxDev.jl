# Trichoplax development model
# MGP Dec 2019

using Makie
using LinearAlgebra
using Statistics

struct Trichoplax
    L0::Float64       # un-stressed cytoskeleton edge length
    V0::Float64       # un-stressed cell volume (area in 2D)
    vertex::Array     # nVertex*2 cell vertices
    skeleton::Array   # nLink*2   indices of 2 cell vertices
    cell::Array       # nCells*6  indices of cell vertices
    skin::Array       # nSkin*1  indices of exterior skeleton links
end

diameter = 36

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
    MaxNCells = Int64(8*ceil((nRows+1)*(nCols+1)))

    # cell centre x-y coordinate array
    cellNucleus = NaN*ones(MaxNCells,2)

    # cell vertex array
    # each row is a vertex (x,y) of a hexagonal cell boundary
    # TODO: Determine empirically how big (or small) this array needs to be
    cellVertex = NaN*ones(6*MaxNCells, 2)

    # cell boundary indices
    # each hexagonal cell defined by 6 vertex indices listed anticlockwise
    iCell= fill(-1, MaxNCells, 6)

    # links
    # each row defines a link between cell vertices
    # by specifying a pair of row indices in the vertex array
    skeleton = fill(-1, 6*MaxNCells, 2)

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
                    iCell[cell, i] = j
                    vertexExists = true
                end
            end
            if !vertexExists
                nCellVertex = nCellVertex + 1
                cellVertex[nCellVertex, :] = newVertex
                iCell[cell, i] = nCellVertex
            end
            # create links
            if i>1
                linkExists = false
                for j in 1:nLink
                    if (transpose([iCell[cell,i-1] iCell[cell,i]]') ==
                                                          skeleton[j,:]') |
                       (transpose([iCell[cell,i] iCell[cell,i-1]]') ==
                                                          skeleton[j,:]')
                        linkExists = true
                    end
                end
                if !linkExists
                    nLink = nLink + 1
                    skeleton[nLink,:] =
                                 [iCell[cell,i-1] iCell[cell,i]]
                end
            end
            if nLink>0
                lines!(cellVertex[skeleton[nLink, :],1],
                    cellVertex[skeleton[nLink, :],2])
                # display(petridish)
                # sleep(1/25)
            end
       end
       # add 6th edge (close the loop from 6th to 1st vertex)
       linkExists = false
       for j in 1:nLink
           if (transpose([iCell[cell,1] iCell[cell,6]]') == skeleton[j,:]') |
              (transpose([iCell[cell,6] iCell[cell,1]]') == skeleton[j,:]')
               linkExists = true
           end
       end
       if !linkExists
           nLink = nLink + 1
           skeleton[nLink,:] = [iCell[cell,1] iCell[cell,6]]
       end
    end

    # find skin segments (they have a vertex with only 2 links)
    # nb skin segments with two exterior vertices are duplicated
    #    duplicates are removed by unique() in Trichoplax constructor call
    nSkin = 0
    skin = fill(0, 4*nCells)
    for iVertex in 1:nCellVertex
        skinCount = 0
        for iLink in 1:nLink  # count links to this vertex
            if any(skeleton[iLink,:].==iVertex)
                nSkin = nSkin + 1
                skinCount = skinCount + 1
                skin[nSkin] = iLink
            end
        end
        if skinCount>2  # not a skin segment
            nSkin = nSkin - skinCount  # reset pointer to overwrite
        end
    end

   # construct Trichoplax
    trichoplax = Trichoplax(cellDiam/2.0,
                            3.0*sqrt(3)*(cellDiam/2.0)^2/2.0,
                            cellVertex[1:nCellVertex,:].*h,
                            skeleton[1:nLink,:],
                            iCell[1:nCells,:],
                            unique(skin[1:nSkin])
                            )

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

    # render cells
    poly!(vertex, facet, color = :lightblue)


    # render skeletons on top
    for link in 1:size(trichoplax.skeleton,1)
        lines!(scene, trichoplax.vertex[trichoplax.skeleton[link, :],1],
               trichoplax.vertex[trichoplax.skeleton[link, :],2])
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
# for i in 1:size(trichoplax.vertex, 1)
#     trichoplax.vertex[i,:] = trichoplax.vertex[i,:] + 0.5*(rand(2).-0.5)
# end




# # enclosing circle
# nPt = 128.
# lines!(diameter/2.0*cos.(2π*(0:nPt)./nPt),
#         diameter/2.0*sin.(2π*(0:nPt)./nPt),
#         color = :lightblue)




# println(MaxNCells)
#
# function cellVolume(cell, trichoplax)
# # volume (area in 2D) of cell
#
#     V = 0.0
#     for i in 1:6
#         j = i % 6 + 1
#         V = V + trichoplax.vertex[trichoplax.cell[cell,i], 1] *
#                 trichoplax.vertex[trichoplax.cell[cell,j], 2] -
#                 trichoplax.vertex[trichoplax.cell[cell,i], 2] *
#                 trichoplax.vertex[trichoplax.cell[cell,j], 1]
#     #    A = A + x[i]*y[j] - y[i]*x[j]
#     end
#
#     return (abs(V/2.0))
#
# end

function cellVolume(x,y)
    # volume (area in 2D) of cell with vertices x (Mx2)

    V = 0.0
    for i in 1:6
        j = i % 6 + 1
        V = V + x[i]*y[j] - x[j]*y[i]
    end

    return (abs(V/2.0))


end


function stress(trichoplax, dx)
    # mechanical stress on trichoplax with each vertex coordinate
    # perturbed in turn  by  dx
    # returned value E is the same size as trichoplax.vertex
    # E[i,1] is the stress when vertex[i,:] (ith x-ccord) is perterbed by dx
    # E[i,2] is the stress when vertex[i,:] (ith y-coord) is perterbed by dx
    #        and all other vertices stay fixed.
    # Used to compute the gradient of stress given the vertices
    # nb this is inefficient because it computes the stress over the whole
    # animal but moving one vertex while holding others fixed affects
    # only neighbouring cells (TBD)

    # copy vertices
    v = copy(trichoplax.vertex)
    E = fill(0.0, size(v)...)
    for i in 1:size(v,1)
        v[i,1] = v[i,1].+dx   # perturb ith vertex in x-direction

        # turgor pressure
        for cell in 1:size(trichoplax.cell,1)
            x = v[trichoplax.cell[cell,:],1]
            y = v[trichoplax.cell[cell,:],2]
            E[i,1] = E[i,1] + (cellVolume(x, y)-trichoplax.V0)^2
        end

        # skeleton tension
        for link in 1:size(trichoplax.skeleton, 1)

            xdiff = v[trichoplax.skeleton[link,1],1] -
                 v[trichoplax.skeleton[link,2],1]
            ydiff = v[trichoplax.skeleton[link,1],2] -
                      v[trichoplax.skeleton[link,2],2]

            E[i,1] = E[i, 1] .+ abs(sqrt(xdiff^2 + ydiff^2) - .9*trichoplax.L0)


        end

        v[i,1] =  v[i,1].-dx  # put the vertex back

        v[i,2] =  v[i,2].+dx   # perturb ith vertex in y-direction

        for cell in 1:size(trichoplax.cell,1)
            x = v[trichoplax.cell[cell,:],1]
            y = v[trichoplax.cell[cell,:],2]
            E[i,2] = E[i,2]+ (cellVolume(x, y)-trichoplax.V0)^2
        end

        # skeleton tension
        for link in 1:size(trichoplax.skeleton,1)

            xdiff = v[trichoplax.skeleton[link,1],1] -
                 v[trichoplax.skeleton[link,2],1]
            ydiff = v[trichoplax.skeleton[link,1],2] -
                 v[trichoplax.skeleton[link,2],2]

            E[i,2] = E[i, 2] .+ abs(sqrt(xdiff^2 + ydiff^2) - .9*trichoplax.L0)

        end

        v[i,2] =  v[i,2].-dx  # put the vertex back

    end

    E
end

for i in 1:size(trichoplax.vertex, 1)
    trichoplax.vertex[i,:] = trichoplax.vertex[i,:] + (rand(2).-.5)
end
draw_trichoplax(trichoplax, petridish)


function morph(trichoplax)
dx = 1e-8
frameCount = 0
for i in 1:1024
    D =  (stress(trichoplax, dx) - stress(trichoplax, -dx))/(2.0*dx)
    #println(D[1,:])
    for j in 1:size(trichoplax.vertex,1)
        #println(D[j,:])
        trichoplax.vertex[j,:] = trichoplax.vertex[j,:]- .001*D[j,:]
    end
    # shift to origin
    x0 = mean(trichoplax.vertex, dims=1)
    #println(x0)
    # for j in 1:size(trichoplax.vertex,1)
    #     trichoplax.vertex[j,:] = trichoplax.vertex[j,:]- transpose(x0)
    # end


    # scene = Scene(limits = FRect(-sceneWidth/2, -sceneWidth/2,
    #         sceneWidth,sceneWidth), scale_plot = false)
    frameCount = frameCount + 1
    if frameCount > 64
    draw_trichoplax(trichoplax, petridish)
  display(petridish)
  sleep(.1)
  frameCount = 0
end
end
end

# dx = 0.1
#
# for i in 1:12
#     for j in 1:size(trichoplax.vertex,1)
#         trichoplax.vertex[j,:] = trichoplax.vertex[j,:]+0.5*(rand(2).-.5)
#     end
#
# scene = Scene(limits = FRect(-sceneWidth/2, -sceneWidth/2,
#         sceneWidth,sceneWidth), scale_plot = false)
# draw_trichoplax(trichoplax, scene)
# display(scene)
# sleep(.1)
# end
