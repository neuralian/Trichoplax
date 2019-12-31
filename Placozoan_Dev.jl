# Placozoan_Dev.jl

using Makie
using Colors
using Placozoan

struct Trichoplax
  layers::Int64
  vertex::Array{Float64}
  cell::Array{Int64}


end

function make_trichoplax(layers, layerwidth, dlayers, mlayers)
  (dvertex,dlink,dlayer) = delaunaydisc(dlayers, layerwidth)
  nbrs=neighbours(dvertex,dlink,dlayer)

  # make celllayers layers of hexagonal cells (Voronoi tesselation)
  ncells = 3*layers*(layers-1)+1
  (vertex, cell, link, layer) = makebody(dvertex, nbrs, dlayer[1:layers+1])

  perimeter = findperimeter(cell, link, layer)
  vertex = smoothperimeter(vertex, perimeter...)
  rfcenter = makereceptivefields(vertex, cell, layer,mlayers);


  (mapvertex,mapcell) = makecellmap(vertex, cell, layer,mlayers);


  return Trichoplax(layers, vertex, cell)
end


# MAIN
dlayers = 8    # number of layers in Delaunay tesselation
layers = 7   # number of body cell layers
mlayers = 2     # map layers
layerwidth = 5.0
(dvertex,dlink,dlayer) = delaunaydisc(dlayers, layerwidth)
nbrs=neighbours(dvertex,dlink,dlayer)

# make celllayers layers of hexagonal cells (Voronoi tesselation)
ncells = 3*layers*(layers-1)+1
(vertex, cell, link, layer) = makebody(dvertex, nbrs, dlayer[1:layers+1])

perimeter = findperimeter(cell, link, layer)
vertex = smoothperimeter(vertex, perimeter...)
rfcenter = makereceptivefields(vertex, cell, layer,mlayers);


(mapvertex,mapcell) = makecellmap(vertex, cell, layer,mlayers);

# Draw
s = Scene(resolution = (800,800), scale_plot = false)

# draw cell nuclei
scatter!(
  dvertex[1:ncells, 1],
  dvertex[1:ncells, 2],
  markersize = layerwidth / 6.0,
  color = RGB(1.0, .5, 0.2),
)

drawcells(vertex,cell, :red, 1.0)
# scatter!(
#   rfcenter[:, 1],
#   rfcenter[:, 2],
#   marker = :circle,
#   strokewidth = layerwidth/6.,
#   strokecolor = :blue,
#   markersize = layerwidth / 4.0,
#   color = :white
# )
drawcells(mapvertex,mapcell, RGB(.25, 0.5,.4))
display(s)
