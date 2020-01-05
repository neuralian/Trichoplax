# Placozoan_Dev.jl

using Makie
using Colors
#using LinearAlgebra
using Statistics
using Placozoan


function make_trichoplax(layers, mapdepth, mat)

  nbrs=neighbours(mat.vertex,mat.edge,mat.layer)

  # make celllayers layers of hexagonal cells (Voronoi tesselation)
  ncells = 3*layers*(layers-1)+1
  (vertex, cell, link, layer) =
                  makebody(mat.vertex, nbrs, mat.layer[1:layers+1])

  perimeter = findperimeter(cell, link, layer)
  vertex = smoothperimeter(vertex, perimeter...)

  (mapvertex,mapcell) = makecellmap(vertex, cell, layer, mapdepth);

  return Trichoplax(layers, vertex, cell, layer, link, perimeter,
                    mapdepth, mapvertex, mapcell)
end


function draw(trichoplax::Trichoplax, color=:black, linewidth = 1.0)
  drawcells(trichoplax.vertex, trichoplax.cell, color, linewidth)
end

function drawmap(trichoplax, color = :blue, linewidth = 0.5)
  drawcells(trichoplax.mapvertex, trichoplax.mapcell, color, linewidth)

end

function draw(tridisc, color = RGB(.5,.5,.5), linewidth = 0.25)
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
worldlayers = 8   # number of layers (rings) in mat
bodylayers = 5   # number of body cell layers
mapdepth = 2     # map layers
layerwidth = 5.0

mat = Tridisc(worldlayers, layerwidth)
trichoplax = make_trichoplax(bodylayers, mapdepth, mat)



# Draw
s = Scene(resolution = (800,800), scale_plot = false)
draw(trichoplax, :red)
drawmap(trichoplax, :orange, 1.0)
draw(mat, RGB(.5, .5, 1.0))

# # draw cell nuclei
# scatter!(
#   dvertex[1:ncells, 1],
#   dvertex[1:ncells, 2],
#   markersize = layerwidth / 6.0,
#   color = RGB(1.0, .5, 0.2),
# )
#
# drawcells(vertex,cell, :red, 1.0)
# # scatter!(
# #   rfcenter[:, 1],
# #   rfcenter[:, 2],
# #   marker = :circle,
# #   strokewidth = layerwidth/6.,
# #   strokecolor = :blue,
# #   markersize = layerwidth / 4.0,
# #   color = :white
# # )
# drawcells(mapvertex,mapcell, RGB(.25, 0.5,.4))
display(s)
