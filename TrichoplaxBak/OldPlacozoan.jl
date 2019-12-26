# Trichoplax development model
# MGP Dec 2019
#module Placozoan

using Makie
using LinearAlgebra
using Statistics
using Colors

#export Trichoplax,  draw, morph


struct Trichoplax
    L0::Float64       # un-stressed cytoskeleton edge length
    V0::Float64         # un-stressed cell volume (area in 2D)
    vertex::Array     # nVertex*2 cell vertices
    skeleton::Array   # nLink*2   indices of 2 cell vertices
    cell::Array       # nCells*6  indices of cell vertices
    skin::Array       # nSkin*1  indices of exterior skeleton links
end

function make_trichoplax(cellDiam, Nlayers, nMorph = 64)
    # Construct Trichoplax with N layers
    # Build hexagonal pre-morph, then morph by minimising
    #  energy = turgor pressure + cytoskeleton elasticity + cell surface energy
    #  nMorph = number of morphing steps

    # cell edge length re. diameter
    h = sqrt(3.0)/2.0

    # vertex distance tolerance
    vertexTol = cellDiam/100.

    maxNCells = Int64(round(4*Nlayers^2 - Nlayers*(Nlayers-1)))

    # x-y coordinates of cell centres
    cellNucleus = fill(0.0, maxNCells, 2)

    # build in 1st quadrant, then copy to other quadrants
    #  (don't copy vertices on axes, ie 1 coordinate is zero)
    nCells = 0
    for i in 1:Nlayers
        y = cellDiam*(i-1)*1.5/2
        for j in 1:Int64(ceil((2Nlayers-i)/2))
            x = h*cellDiam*(j - 0.5 - 0.5*isodd(i))
            nCells = nCells+1
            cellNucleus[nCells, :] = [x y]
            if i>1
                nCells = nCells + 1
                cellNucleus[nCells,:] = [x -y]
            end
            if x>cellDiam/4
                nCells = nCells + 1
                cellNucleus[nCells, :] = [-x y]
            end
            if (x > cellDiam/4.) &  (y > cellDiam/4.)
                nCells = nCells + 1
                cellNucleus[nCells,:] = [-x -y]
            end
        end
    end

    # draw nuclei - for debug; comment out to run
    scatter!(cellNucleus[:,1], cellNucleus[:,2], markersize = cellDiam/5)

    # cell vertex array
    # each row is a vertex (x,y) of a hexagonal cell boundary
    # TODO: Determine empirically how big (or small) this array needs to be
    cellVertex = NaN*ones(6*nCells, 2)

    # cell boundary indices
    # each hexagonal cell defined by 6 vertex indices listed anticlockwise
    iCell= fill(-1, nCells, 6)

    # links
    # each row defines a link between cell vertices
    # by specifying a pair of row indices in the vertex array
    skeleton = fill(-1, 6*nCells, 2)

    # construct cells and links
    nCellVertex = 0
    nLink = 0
    for cell in 1:nCells
        for i in 1:6
            newVertex = cellNucleus[cell,:]' .+
                         0.5*cellDiam.*[cos(2π*(i-.5)/6)  sin(2π*(i-.5)/6)]
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

            # # draw links - for debug; comment out to run
            # if nLink>0
            #     lines!(cellVertex[skeleton[nLink, :],1],
            #         cellVertex[skeleton[nLink, :],2])
            # end

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
    skin = fill(0, 8*nCells)
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


   V0 = 3.0*sqrt(3)*(cellDiam/2.0)^2/2.0


   # construct Trichoplax
    trichoplax = Trichoplax(cellDiam/2.0,
                            V0,
                            cellVertex[1:nCellVertex,:],
                            skeleton[1:nLink,:],
                            iCell[1:nCells,:],
                            unique(skin[1:nSkin])
                            )


    morph(trichoplax, nMorph)


    # return
    trichoplax

end



function draw(trichoplax)

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
    poly!(vertex, facet, color = RGBA(0.75, .25, .25, .5))


    # render skeletons on top
    for link in 1:size(trichoplax.skeleton,1)
        lines!(trichoplax.vertex[trichoplax.skeleton[link, :],1],
               trichoplax.vertex[trichoplax.skeleton[link, :],2])
    end

end

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

    c = 1.0  # skeleton pre-tension
    σ = 100.0  # surface energy density

    # copy vertices
    v = copy(trichoplax.vertex)
    E = fill(0.0, size(v)...)
    for i in 1:size(v,1)
        v[i,1] = v[i,1].+dx   # perturb ith vertex in x-direction

        # cell turgor pressure
        for cell in 1:size(trichoplax.cell,1)
            x = v[trichoplax.cell[cell,:],1]
            y = v[trichoplax.cell[cell,:],2]
#            E[i,1] = E[i,1] + log(cellVolume(x, y)/trichoplax.V0)^2
            E[i,1] = E[i,1] + (cellVolume(x, y)- trichoplax.V0)^2
        end

        # cytoskeleton tension
        for link in 1:size(trichoplax.skeleton, 1)
            xdiff = v[trichoplax.skeleton[link,1],1] -
                 v[trichoplax.skeleton[link,2],1]
            ydiff = v[trichoplax.skeleton[link,1],2] -
                      v[trichoplax.skeleton[link,2],2]
            E[i,1] = E[i, 1] .+ log(sqrt(xdiff^2 + ydiff^2)/trichoplax.L0*c)^2
        end

        # skin surface tension
        for link in 1:size(trichoplax.skin,1)
            xdiff = v[trichoplax.skeleton[trichoplax.skin[link],1],1] -
                    v[trichoplax.skeleton[trichoplax.skin[link],2],1]
            ydiff = v[trichoplax.skeleton[trichoplax.skin[link],1],2] -
                    v[trichoplax.skeleton[trichoplax.skin[link],2],2]
            E[i,1] = E[i,1] .+ σ*sqrt(xdiff^2 + ydiff^2)
        end

        v[i,1] =  v[i,1].-dx  # put the x-vertex back

        v[i,2] =  v[i,2].+dx   # perturb ith vertex in y-direction

        # cell turgor pressure
        for cell in 1:size(trichoplax.cell,1)
            x = v[trichoplax.cell[cell,:],1]
            y = v[trichoplax.cell[cell,:],2]
            E[i,2] = E[i,2]+ (cellVolume(x, y)-trichoplax.V0)^2
        end

        # cytoskeleton tension
        for link in 1:size(trichoplax.skeleton,1)
            xdiff = v[trichoplax.skeleton[link,1],1] -
                 v[trichoplax.skeleton[link,2],1]
            ydiff = v[trichoplax.skeleton[link,1],2] -
                 v[trichoplax.skeleton[link,2],2]
            E[i,2] = E[i, 2] .+ log(sqrt(xdiff^2 + ydiff^2)/trichoplax.L0*c)^2

        end

        # skin surface tension
        for link in 1:size(trichoplax.skin,1)
            xdiff = v[trichoplax.skeleton[trichoplax.skin[link],1],1] -
                    v[trichoplax.skeleton[trichoplax.skin[link],2],1]
            ydiff = v[trichoplax.skeleton[trichoplax.skin[link],1],2] -
                    v[trichoplax.skeleton[trichoplax.skin[link],2],2]
            E[i,2] = E[i,2] .+ σ*sqrt(xdiff^2 + ydiff^2)
        end

        v[i,2] =  v[i,2].-dx  # put the vertex back

    end

    E
end



function morph(trichoplax, N=16, dx = 1e-2)
    # morph to minimise stress
    frameCount = 0
    for i in 1:N
        D =  (stress(trichoplax, dx) - stress(trichoplax, -dx))/(2.0*dx)
        #println(D[1,:])
        for j in 1:size(trichoplax.vertex,1)
            #println(D[j,:])
            trichoplax.vertex[j,:] = trichoplax.vertex[j,:]- .01*D[j,:]
        end
        # shift to origin
        x0 = mean(trichoplax.vertex, dims=1)
       println(i)
    end
end


#end  # module Placozoan


Nlayers =4
cellDiam = 7.5


# Set scene
sceneWidth = Int64(round(cellDiam*Nlayers*1.5))
if isodd(sceneWidth)
    sceneWidth = sceneWidth + 1
end

world = Scene(limits = FRect(-sceneWidth/2, -sceneWidth/2,
                              sceneWidth,sceneWidth), scale_plot = false)

@time trichoplax = make_trichoplax(cellDiam, Nlayers, 32);
draw(trichoplax)
display(world)
println(size(trichoplax.skin))
