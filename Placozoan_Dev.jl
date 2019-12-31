# Placozoan_Dev.jl

using Makie
using Colors
using Placozoan

struct Mat
  layers::Int64
  vertex::Array{Float64}
  link::Array{Int64}
end

struct Trichoplax
  layers::Int64                      # number of layers of cells
  vertex::Array{Float64}             # cell vertices
  cell::Array{Int64}                 # [i,j] index jth vertex of ith cell
  layer::Vector{Array{Int64,1} }     # [i][j] index jth cell in ith layer
  link::Array{Int64}                 # [i,:] index links between cells
  perimeter                          #  perimeter vertex indices (3-tuple)
  mapdepth::Int64                    # number of cell layers in sensory map
end

function make_trichoplax(layers, layerwidth, dlayers, mapdepth)

  (dvertex,dlink,dlayer) = delaunaydisc(dlayers, layerwidth)
  nbrs=neighbours(dvertex,dlink,dlayer)

  # make celllayers layers of hexagonal cells (Voronoi tesselation)
  ncells = 3*layers*(layers-1)+1
  (vertex, cell, link, layer) = makebody(dvertex, nbrs, dlayer[1:layers+1])

  perimeter = findperimeter(cell, link, layer)
  vertex = smoothperimeter(vertex, perimeter...)
  #rfcenter = makereceptivefields(vertex, cell, layer,mapdepth);


  (mapvertex,mapcell) = makecellmap(vertex, cell, layer,mapdepth);


  return Trichoplax(layers, vertex, cell, layer, link, perimeter, mapdepth)
end


# MAIN
dlayers = 5    # number of layers in Delaunay tesselation
layers = 4   # number of body cell layers
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

trichoplax = make_trichoplax(layers, layerwidth, dlayers, mapdepth)

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
