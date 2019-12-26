# RoundDelaunay.jl
# circular Delaunay triangulation

using Makie


<<<<<<< Updated upstream
function delaunayDisc(nLayers, layerWidth)
    # returns (cellNucleus, edge, layer)
    # cellnucleus = nCells x 2 coordinates of ith cell nucleus
    # edge = nEdge x indices of cell nuclei defining the edge
    # layer = ith row indexes cell nuclei in each layer
    # MGP Dec 2019
=======
nLayers = 12
layerW = 5.0
>>>>>>> Stashed changes
nCells = Int64(round(π*(nLayers+1)^2))
cellNucleus = fill(0.0,nCells , 2)
edge = fill(0, 6*nCells, 2)
Layer = fill(Int64[],nLayers)
nLayer = vcat([1], [6*i for i in 1:nLayers-1]) # n cells in layer
nCLink = fill(0, nCells, 1)  # number of Delaunay links per cell
iCell = 1
cellNucleus[iCell,:] = [0.0 0.0]  # cell at origin
Layer[1,1] = [1]
iEdge =0
w = 0
for iLayer in 2:nLayers
    Layer[iLayer,:] = [fill(0,nLayer[iLayer])]
    w = w + 2π*(rand(1)[1] - 0.5)/nLayer[iLayer]
    for j in 1:nLayer[iLayer]
        iCell = iCell + 1
        Layer[iLayer][j] = iCell
        cellNucleus[iCell,:] =
        (iLayer-1)*layerWidth.*[cos(2π*((j-1)/nLayer[iLayer])+w)
                       sin(2π*((j-1)/nLayer[iLayer])+w)]

        # Delaunay triangulation
        if j>1
            iEdge = iEdge + 1
            edge[iEdge, :] = [iCell-1 iCell]
            nCLink[iCell-1] = nCLink[iCell-1]+1
            nCLink[iCell] = nCLink[iCell]+1
            #lines!(cellNucleus[edge[iEdge,:],1], cellNucleus[edge[iEdge,:],2])
            #display(s)
        end
    end
    iEdge = iEdge + 1
    edge[iEdge, :] = [iCell Layer[iLayer][1]]
    nCLink[iCell] = nCLink[iCell]+1
    nCLink[Layer[iLayer][1]] = nCLink[Layer[iLayer][1]]+1
    #lines!(cellNucleus[edge[iEdge,:],1], cellNucleus[edge[iEdge,:],2])
    #display(s)
end
nCells = iCell
cellNucleus = cellNucleus[1:nCells,:]

# scatter!(cellNucleus[:,1], cellNucleus[:,2], markersize = layerWidth/4.)
# for i in 1:nLayers
#      scatter!([cellNucleus[Layer[i][1], 1]], [cellNucleus[Layer[i][1],2]],
#      markersize = .5, color = :red);
#  end
# display(s)

# complete Delaunay triangles
for j in 1:6      # center cell special case, 6 links
    iEdge = iEdge + 1
    edge[iEdge,:] = [1 j+1]
    nCLink[1] = nCLink[1]+1
    nCLink[j+1] = nCLink[j+1]+1
    #lines!(cellNucleus[edge[iEdge,:],1], cellNucleus[edge[iEdge,:],2])
    #display(s)
end

iLink = []
Link = []
for iLayer in 2:nLayers-1
    iLink = [ nLayer[iLayer+1] 1 2]
    Link = Layer[iLayer+1][iLink[1:3]]
    for i in 1:length(Layer[iLayer])
        #println((iLayer, Link))
        for j in 1:length(Link)            #global iEdge = iEdge + 1
            iEdge = iEdge + 1
            edge[iEdge, :] = [Layer[iLayer][i] Link[j]]
            nCLink[Layer[iLayer][i]] = nCLink[Layer[iLayer][i]]+1
            nCLink[Link[j]] = nCLink[Link[j]]+1
            # lines!(cellNucleus[edge[iEdge,:],1], cellNucleus[edge[iEdge,:],2])
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
edge = edge[1:iEdge,:]

return (cellNucleus, edge, Layer)

end


s = Scene(scale_plot = false)
nLayer = 4
layerWidth = 5.0
DD = delaunayDisc(nLayer, layerWidth)
cellNucleus=DD[1]; edge = DD[2];  layer = DD[3];

# cell nucleuslei
scatter!(cellNucleus[:,1], cellNucleus[:,2], markersize = layerWidth/4.)
# Delaunay
for i in 1:size(edge,1)
  lines!(cellNucleus[edge[i,:],1], cellNucleus[edge[i,:],2])
end
for i in 1:nLayer
     scatter!([cellNucleus[layer[i][1], 1]], [cellNucleus[layer[i][1],2]],
     markersize = .5, color = :red);
 end
display(s)
