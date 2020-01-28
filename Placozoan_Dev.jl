# Placozoan_Dev.jl

using Makie
using Colors
#using LinearAlgebra
using Statistics
using Placozoan


function draw(trichoplax::Trichoplax, color=:black, linewidth = 1.0)
  drawcells(trichoplax.vertex, trichoplax.cell, color, linewidth)
end

function drawmap(trichoplax, color = :blue, linewidth = 0.5)
  drawcells(trichoplax.mapvertex, trichoplax.mapcell, color, linewidth)

end

function draw(tridisc::Tridisc, color = RGB(.5,.5,.5), linewidth = 0.25)
  for link in 1:size(tridisc.edge,1)
      lines!(tridisc.vertex[tridisc.edge[link, :],1],
             tridisc.vertex[tridisc.edge[link, :],2],
             color=color, linewidth=linewidth)
  end
end


function paint(trichoplax, color = :red)

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
    poly!(vertex, facet, color = rand(size(vertex,1)),
            colormap = :blues)

    #
    # # render skeletons on top
    # for link in 1:size(trichoplax.skeleton,1)
    #     lines!(trichoplax.vertex[trichoplax.skeleton[link, :],1],
    #            trichoplax.vertex[trichoplax.skeleton[link, :],2])
    # end

end

# MAIN
bodylayers = 12   # number of body cell layers
mapdepth = 4     # map layers
celldiam = 10.0


trichoplax = Trichoplax(bodylayers, celldiam, mapdepth)

# Draw
s = Scene(resolution = (800,800), scale_plot = false)
draw(trichoplax, :red)
#drawmap(trichoplax, :orange, 1.0)
draw(trichoplax.skeleton, RGB(.5, .5, 1.0))

display(s)
