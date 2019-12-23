# Placozoan_Dev.jl

using Makie
using Colors
using Placozoan

function neighbours(nucleus, dlink, layer)
  # index the six neighbour nuclei of each nucleus (excepting outer layer)
  # given output from delaunayDisc()...   (nb splatted argument)
  # outer nucleus layer
  # MGP Dec 2019

  nCells = size(nucleus, 1)
  neighbour = fill(0, nCells, 6)   # neighbour indices for each cell
  nbrCount = fill(0, nCells)       # number of neighbours found

  # find neighbours of each cell
  for i in 1:size(dlink,1)
    a = dlink[i,1]
    b = dlink[i,2]
    nbrCount[a] = nbrCount[a]+1
    nbrCount[b] = nbrCount[b]+1
    neighbour[a, nbrCount[a]] = b
    neighbour[b, nbrCount[b]] = a
  end
  # put neighbours in counterclockwise order
  # nb row 1 is already ccw, skip external nuclei (less than 6 links)
  for i in 2:(nCells-length(layer[end]))
    dx = nucleus[neighbour[i,:],1].-nucleus[i,1]
    dy = nucleus[neighbour[i,:],2].-nucleus[i,2]
    order = sortperm(atan.(dy,dx))
    neighbour[i,:] = neighbour[i, order]
  end
  neighbour
end

function cells(nucleus, neighbour, layer)
  # construct cell vertices
  # MGP Dec 2019

  nCells = size(nucleus,1) - length(layer[end])
  vertex = fill(0.0, nCells, 6, 2 )  # 6 vertices per cell

  for i in 1:nCells
    x0 = nucleus[i,1]
    y0 = nucleus[i,2]
    for j in 1:6
        k = neighbour[i,j]
        m = neighbour[i, j%6+1]
        vertex[i,j, :] = [(x0 + nucleus[k,1] + nucleus[m,1])/3.,
                          (y0 + nucleus[k,2]+ nucleus[m,2])/3.]
    end
  end
  vertex
end

# function drawcell(iCell, cell, vertex)
#  lines!(vtx[iCell,[1:end; 1],1], vtx[iCell,[1:end; 1],2])
# end

nLayer = 8
layerWidth = 5.0
@time DD = delaunayDisc(nLayer, layerWidth)
@time nbrs=neighbours(DD...)


cellNucleus = DD[1];
dlink = DD[2];
layer = DD[3];

@time vtx = cells(cellNucleus, nbrs, layer)

# Draw
s = Scene(resolution = (800,800), scale_plot = false)
# cell nucleuslei
scatter!(
  cellNucleus[:, 1],
  cellNucleus[:, 2],
  markersize = layerWidth / 6.0,
  color = RGB(0.7, 0.7, 0.7),
)
# Delaunay
@inbounds for i = 1:size(dlink, 1)
  lines!(
    cellNucleus[dlink[i, :], 1],
    cellNucleus[dlink[i, :], 2],
    color = RGB(0.7, 0.7, 0.7), linewidth = 0.25
  )
end
@inbounds for i = 1:nLayer
  scatter!(
    [cellNucleus[layer[i][1], 1]],
    [cellNucleus[layer[i][1], 2]],
    markersize = layerWidth / 12,
    color = :red,
  )
end


for i in 1:size(vtx,1)
  lines!(vtx[i,[1:end; 1],1], vtx[i,[1:end; 1],2])
end
display(s)
