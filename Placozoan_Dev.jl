# Placozoan_Dev.jl

using Makie
using Colors
using Placozoan

function neighbours(nucleus, dlink, layer)
  # index the six neighbour nuclei of each nucleus (excepting outer layer)
  # given output from delaunayDisc()...   (nb splatted argument)
  # outer nucleus layer
  # MGP Dec 2019

  nLayer = size(layer,1)
  N = [size(layer[i], 1) for i in 1:nLayer]         # nuclei per layer
  nCells = size(nucleus, 1) - length(layer[nLayer]) # 1 less layer of cells
  neighbour = fill(0, nCells, 6)
  nD = size(dlink,1)  # number of delaunay links

  # Centre cell (layer 1) is a special case
  # its neighbours are the 6 cells in layer 2
  neighbour[1,:] = collect(2:7)

  for i = 2:(nLayer-1)       # for each layer
    iCandy = -1              # offset to next candidate neighbour in next layer
    for j = 1:N[i]           # for each cell in layer
      # first neighbour is previous cell in layer
      neighbour[layer[i][j], 1] = layer[i][mod(j-2,N[i])+1]
      # next candidate neighbour is previous cell in next layer
      candidate = layer[i+1][mod(j+iCandy-1,N[i+1])+1]
      println((i,j, candidate))
      nNb = 0   # number of neighbours found in outer layer
      for k in 1:3  # there are 2 or 3 neighbours in the next layer
        for id = 1:nD
          # if there is a d-link from current cell to candidate
          if any(dlink[id,:].==layer[i][j]) & any(dlink[id,:].==candidate)
            nNb = nNb + 1
            neighbour[layer[i][j], 1 + nNb] = candidate
            candidate = layer[i+1][mod(j+iCandy-1,N[i+1])+1] # next candidate
          end
        end
      end
      iCandy = iCandy + nNb
    end
  end
  neighbour
end




nLayer = 4
layerWidth = 5.0
@time DD = delaunayDisc(nLayer, layerWidth)
cellNucleus = DD[1];
dlink = DD[2];
layer = DD[3];

nb=neighbours(DD...)
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
    color = RGB(0.7, 0.7, 0.7),
  )
end
@inbounds for i = 1:nLayer
  scatter!(
    [cellNucleus[layer[i][1], 1]],
    [cellNucleus[layer[i][1], 2]],
    markersize = layerWidth / 6,
    color = :red,
  )
end
display(s)
