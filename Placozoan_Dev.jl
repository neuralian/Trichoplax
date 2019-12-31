# Placozoan_Dev.jl

using Makie
using Colors
using Placozoan



# MAIN
dlayers = 8    # number of layers in Delaunay tesselation
blayers = 6   # number of body cell layers
maplayers = 2     # map layers
layerwidth = 5.0
(dvertices,dlink,dlayer) = delaunaydisc(dlayers, layerwidth)
nbrs=neighbours(dvertices,dlink,dlayer)

# make celllayers layers of hexagonal cells (Voronoi tesselation)
ncells = 3*blayers*(blayers-1)+1
(vertex, cell, vlink, clayer) = makebody(dvertices, nbrs, dlayer[1:blayers+1])

perimeter = findperimeter(cell, vlink, clayer)
vertex = smoothperimeter(vertex, perimeter...)
rfcenter = makereceptivefields(vertex, cell, clayer,maplayers);


(cellmapvertex,cellmapcell) = makecellmap(vertex, cell, clayer,maplayers);

# Draw
s = Scene(resolution = (800,800), scale_plot = false)

# draw cell nuclei
scatter!(
  dvertices[1:ncells, 1],
  dvertices[1:ncells, 2],
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
drawcells(cellmapvertex,cellmapcell, RGB(.25, 0.5,.4))
display(s)
