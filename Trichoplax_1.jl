using Placozoan
using Makie

# for make_trichoplax
using LinearAlgebra

# Placozoan package must be on load path.
# In Atom:
#    right-click on folder, select "work in folder"
#    push!(LOAD_PATH, pwd())

function make_trichoplax(Nlayers)
    # Construct Trichoplax with N layers
    # Build hexagonal pre-morph, then morph by minimising
    #  energy = turgor pressure + cytoskeleton elasticity + cell surface energy

    cellDiam = 5.0   # cell diameter microns

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
        y = cellDiam*(i-1)*(1+h)/2.
        for j in 1:Int64(ceil((2Nlayers-i)/2))
            x = cellDiam*(j - 0.5 - 0.5*isodd(i))
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

    # cell nuclei
    # used to construct cell membranes/cells, not members of Trichoplax
    # nCells = 0
    # for row in 0:nRows
    #     for col in 0:nCols
    #         if iseven(row)
    #             x = col*cellDiam
    #             y = h*row*cellDiam
    #             if (x^2 + y^2) <= (plaxDiam/2.)^2
    #                 nCells = nCells + 1
    #                 cellNucleus[nCells,:] = [x y]
    #                 if y>cellDiam/4.
    #                     nCells = nCells + 1
    #                     cellNucleus[nCells,:] = [x -y]
    #                 end
    #                 if x > cellDiam/4.
    #                     nCells = nCells + 1
    #                     cellNucleus[nCells,:] = [-x y]
    #                 end
    #                     if (x > cellDiam/4.) &  (y > cellDiam/4.)
    #                     nCells = nCells + 1
    #                     cellNucleus[nCells,:] = [-x -y]
    #                 end
    #             end
    #         else
    #             x = (col + 0.5)*cellDiam
    #             y = h*row*cellDiam
    #             if (x^2 + y^2) <= (plaxDiam/2.)^2
    #                 nCells = nCells + 1
    #                 cellNucleus[nCells,:] = [x y]
    #                 nCells = nCells + 1
    #                 cellNucleus[nCells,:] = [-x -y]
    #                 nCells = nCells + 1
    #                 cellNucleus[nCells,:] = [x -y]
    #                 nCells = nCells + 1
    #                 cellNucleus[nCells,:] = [-x y]
    #             end
    #         end
    #     end
    # end





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
                            cellVertex[1:nCellVertex,:].*h,
                            skeleton[1:nLink,:],
                            iCell[1:nCells,:],
                            unique(skin[1:nSkin])
                            )

  morph(trichoplax)

  # return
  trichoplax

end


Nlayers = 3
trichoplax = make_trichoplax(Nlayers);


# Set scene
sceneWidth = Int64(round(diameter*1.5))
if isodd(sceneWidth)
    sceneWidth = sceneWidth + 1
end

world = Scene(limits = FRect(-sceneWidth/2, -sceneWidth/2,
                              sceneWidth,sceneWidth), scale_plot = false)


draw(trichoplax)
display(world)
