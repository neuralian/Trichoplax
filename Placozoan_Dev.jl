# Placozoan_Dev.jl

using Makie
using Colors
#using Placozoan

using Distributions
Ndist = Normal()

function delaunayDisc(nLayers, layerWidth)
    # returns (nucleus, dlink, layer)
    # nucleus = nCells x 2 coordinates of ith cell nucleus
    # dlink = ndlink x indices of cell nuclei defining the dlink
    # layer = ith row indexes cell nuclei in each layer
    # MGP Dec 2019
a = .0025   # noise on cell location (re. layerWidth)
nCells = Int64(round(π*(nLayers+1)^2))
nucleus = fill(0.0,nCells , 2)
dlink = fill(0, 6*nCells, 2)
Layer = fill(Int64[],nLayers)
nLayer = vcat([1], [6*i for i in 1:nLayers-1]) # n cells in layer
nCLink = fill(0, nCells, 1)  # number of Delaunay links per cell
iCell = 1
nucleus[iCell,:] = [0.0 0.0]  # cell at origin
Layer[1,1] = [1]
idlink =0
w = 0
@inbounds for iLayer in 2:nLayers
    Layer[iLayer,:] = [fill(0,nLayer[iLayer])]
    w = w + 2.0*(rand(1)[1] - 0.5)/nLayer[iLayer]
    @inbounds for j in 1:nLayer[iLayer]
        dw = rand(Ndist)
        while abs(dw)>3.0
            dw = rand(Ndist)
        end
        dw= .1*dw/(iLayer*layerWidth)
        dr = rand(Ndist)
        while abs(dr)>3.0
            dr = rand(Ndist)
        end
        dr = a*dr
        iCell = iCell + 1
        Layer[iLayer][j] = iCell
        nucleus[iCell,:] =
        (iLayer-1)*(1.0+dr)*layerWidth.*
        [ cos(2π*((j-1)/nLayer[iLayer])+w+dw)
          sin(2π*((j-1)/nLayer[iLayer])+w+dw) ]

        # Delaunay triangulation
        if j>1
            idlink = idlink + 1
            dlink[idlink, :] = [iCell-1 iCell]
            nCLink[iCell-1] = nCLink[iCell-1]+1
            nCLink[iCell] = nCLink[iCell]+1
            #lines!(nucleus[dlink[idlink,:],1], nucleus[dlink[idlink,:],2])
            #display(s)
        end
    end
    idlink = idlink + 1
    dlink[idlink, :] = [iCell Layer[iLayer][1]]
    nCLink[iCell] = nCLink[iCell]+1
    nCLink[Layer[iLayer][1]] = nCLink[Layer[iLayer][1]]+1
    #lines!(nucleus[dlink[idlink,:],1], nucleus[dlink[idlink,:],2])
    #display(s)
end
nCells = iCell
nucleus = nucleus[1:nCells,:]

# scatter!(nucleus[:,1], nucleus[:,2], markersize = layerWidth/4.)
# @inbounds for i in 1:nLayers
#      scatter!([nucleus[Layer[i][1], 1]], [nucleus[Layer[i][1],2]],
#      markersize = .5, color = :red);
#  end
# display(s)

# complete Delaunay triangles
@inbounds for j in 1:6      # center cell special case, 6 links
    idlink = idlink + 1
    dlink[idlink,:] = [1 j+1]
    nCLink[1] = nCLink[1]+1
    nCLink[j+1] = nCLink[j+1]+1
    #lines!(nucleus[dlink[idlink,:],1], nucleus[dlink[idlink,:],2])
    #display(s)
end

iLink = []
Link = []
@inbounds for iLayer in 2:nLayers-1
    iLink = [ nLayer[iLayer+1] 1 2]
    Link = Layer[iLayer+1][iLink[1:3]]
    @inbounds for i in 1:length(Layer[iLayer])
        #println((iLayer, Link))
        @inbounds for j in 1:length(Link)            #global idlink = idlink + 1
            idlink = idlink + 1
            dlink[idlink, :] = [Layer[iLayer][i] Link[j]]
            nCLink[Layer[iLayer][i]] = nCLink[Layer[iLayer][i]]+1
            nCLink[Link[j]] = nCLink[Link[j]]+1
            # lines!(nucleus[dlink[idlink,:],1], nucleus[dlink[idlink,:],2])
            # display(s)
            # sleep(.05)
        end
        #println((Link[end],nCLink[Layer[iLayer][i+1]] ))
        if i<length(Layer[iLayer])
            Link = Link[end] .+ collect(0:(5-nCLink[Layer[iLayer][i+1]]))
        end
        #println(Link)
    end
end
dlink = dlink[1:idlink,:]

return (nucleus, dlink, Layer)

end


function neighbours(nucleus, dlink, layer)
  # index the six neighbour nuclei of each nucleus (excepting outer layer)
  # given output from delaunayDisc()...   (nb splatted argument)
  # outer nucleus layer
  # MGP Dec 2019

  nCells = size(nucleus, 1)
  neighbour = fill(0, nCells, 6)   # neighbour indices for each cell
  nbrCount = fill(0, nCells)       # number of neighbours found

  # find neighbours of each cell
  a = 0
  b = 0
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

  @inbounds for i in 1:nCells
    # x0 = nucleus[i,1]
    # y0 = nucleus[i,2]
     @inbounds for j in 1:6
        k = neighbour[i,j]
        m = neighbour[i, j%6+1]
        vertex[i,j, :] = [(nucleus[i,1] + nucleus[k,1] + nucleus[m,1])/3.,
                          (nucleus[i,2] + nucleus[k,2]+ nucleus[m,2])/3.]

    end
  end
  vertex
end

function drawDelaunayDisc(dlink)
    @inbounds for i = 1:size(dlink, 1)
      lines!(
        cellNucleus[dlink[i, :], 1],
        cellNucleus[dlink[i, :], 2],
        color = RGB(0.7, 0.7, 0.7), linewidth = 0.25
        )
    end
end

# function drawcell(iCell, cell, vertex)
#  lines!(vtx[iCell,[1:end; 1],1], vtx[iCell,[1:end; 1],2])
# end

nLayer = 4
layerWidth = 5.0
print("DelaunayDisc:")
@time DD = delaunayDisc(nLayer, layerWidth)
print("neighbours:")
@time nbrs=neighbours(DD...)


cellNucleus = DD[1];
dlink = DD[2];
layer = DD[3];
print("cells:")
@time vtx = cells(cellNucleus, nbrs, layer)

# Draw
s = Scene(resolution = (800,800), scale_plot = false)
# cell nucleuslei
scatter!(
  cellNucleus[1:(size(cellNucleus,1)-size(layer[nLayer],1)), 1],
  cellNucleus[1:(size(cellNucleus,1)-size(layer[nLayer],1)), 2],
  markersize = layerWidth / 6.0,
  color = RGB(0.7, 0.7, 0.7),
)
# Delaunay

# @inbounds for i = 1:(nLayer-1
#   scatter!(
#     [cellNucleus[layer[i][1], 1]],
#     [cellNucleus[layer[i][1], 2]],
#     markersize = layerWidth / 12,
#     color = :red,
#   )
# end


for i in 1:size(vtx,1)
  lines!(vtx[i,[1:end; 1],1], vtx[i,[1:end; 1],2])
end
display(s)
