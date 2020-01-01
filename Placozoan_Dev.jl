# Placozoan_Dev.jl

using Makie
using Colors
using Placozoan



function make_trichoplax(layers, mapdepth, mat)

#  (dvertex,dlink,dlayer) = delaunaydisc(dlayers, layerwidth)
  nbrs=neighbours(mat.vertex,mat.edge,mat.layer)

  # make celllayers layers of hexagonal cells (Voronoi tesselation)
  ncells = 3*layers*(layers-1)+1
  (vertex, cell, link, layer) =
                  makebody(mat.vertex, nbrs, mat.layer[1:layers+1])

  perimeter = findperimeter(cell, link, layer)
  vertex = smoothperimeter(vertex, perimeter...)
  #rfcenter = makereceptivefields(vertex, cell, layer,mapdepth);


  (mapvertex,mapcell) = makecellmap(vertex, cell, layer, mapdepth);


  return Trichoplax(layers, vertex, cell, layer, link, perimeter, mapdepth)
end


# MAIN
worldlayers = 5    # number of layers (rings) in mat
bodylayers = 4   # number of body cell layers
mapdepth = 1     # map layers
layerwidth = 5.0
# (dvertex,dlink,dlayer) = delaunaydisc(dlayers, layerwidth)
# nbrs=neighbours(dvertex,dlink,dlayer)
#
# # make celllayers layers of hexagonal cells (Voronoi tesselation)
# ncells = 3*layers*(layers-1)+1
# (vertex, cell, link, layer) = makebody(dvertex, nbrs, dlayer[1:layers+1])
#
# perimeter = findperimeter(cell, link, layer)
# vertex = smoothperimeter(vertex, perimeter...)
# rfcenter = makereceptivefields(vertex, cell, layer,mlayers);
#
#
# (mapvertex,mapcell) = makecellmap(vertex, cell, layer,mlayers);


mat = discworld(worldlayers, layerwidth)
trichoplax = make_trichoplax(bodylayers, mapdepth, mat)

function draw(trichoplax::Trichoplax, color=:black)
  drawcells(trichoplax.vertex, trichoplax.cell, color)
end

# Draw
s = Scene(resolution = (800,800), scale_plot = false)
draw(trichoplax, :red)

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
