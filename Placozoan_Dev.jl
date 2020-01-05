# Placozoan_Dev.jl

using Makie
using Colors
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


# MAIN
worldlayers = 16    # number of layers (rings) in mat
bodylayers = 8   # number of body cell layers
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
