using Makie
s = Scene(scale_plot = false)


nLayers = 32
layerW = 5.0
nCells = Int64(round(π*(nLayers+1)^2))
nuc = fill(0.0,nCells , 2)
Dlny = fill(0, 6*nCells, 2)
Layer = fill(Int64[],nLayers+1)
nLayer = vcat([1], [6*i for i in 1:nLayers-1]) # n cells in layer
nCLink = fill(0, nCells, 1)  # number of Delaunay links per cell
iCell = 1
nuc[iCell,:] = [0.0 0.0]  # cell at origin
Layer[1,1] = [1]
iDlny =0
w = 0
for iLayer in 2:nLayers
    global iDlny
    #nLayer = Int64(round(2π*iLayer))
    Layer[iLayer,:] = [fill(0,nLayer[iLayer])]
    global w = w + 2π*(rand(1)[1] - 0.5)/nLayer[iLayer]
    for j in 1:nLayer[iLayer]
        global iCell = iCell + 1
        Layer[iLayer][j] = iCell
        nuc[iCell,:] =
        (iLayer-1)*layerW.*[cos(2π*((j-1)/nLayer[iLayer])+w)
                       sin(2π*((j-1)/nLayer[iLayer])+w)]

        # Delaunay triangulation
        if j>1
            iDlny = iDlny + 1
            Dlny[iDlny, :] = [iCell-1 iCell]
            nCLink[iCell-1] = nCLink[iCell-1]+1
            nCLink[iCell] = nCLink[iCell]+1
            #lines!(nuc[Dlny[iDlny,:],1], nuc[Dlny[iDlny,:],2])
            #display(s)
        end
    end
    iDlny = iDlny + 1
    Dlny[iDlny, :] = [iCell Layer[iLayer][1]]
    nCLink[iCell] = nCLink[iCell]+1
    nCLink[Layer[iLayer][1]] = nCLink[Layer[iLayer][1]]+1
    #lines!(nuc[Dlny[iDlny,:],1], nuc[Dlny[iDlny,:],2])
    #display(s)
end
nCells = iCell
nuc = nuc[1:nCells,:]

# scatter!(nuc[:,1], nuc[:,2], markersize = layerW/4.)
# for i in 1:nLayers
#      scatter!([nuc[Layer[i][1], 1]], [nuc[Layer[i][1],2]],
#      markersize = .5, color = :red);
#  end
# display(s)

# complete Delaunay triangles
for j in 1:6      # center cell special case, 6 links
    global iDlny = iDlny + 1
    Dlny[iDlny,:] = [1 j+1]
    nCLink[1] = nCLink[1]+1
    nCLink[j+1] = nCLink[j+1]+1
    #lines!(nuc[Dlny[iDlny,:],1], nuc[Dlny[iDlny,:],2])
    #display(s)
end

iLink = []
Link = []
for iLayer in 2:nLayers-1
    global iLink = [ nLayer[iLayer+1] 1 2]
    global Link = Layer[iLayer+1][iLink[1:3]]
    for i in 1:length(Layer[iLayer])
        #println((iLayer, Link))
        for j in 1:length(Link)            #global iDlny = iDlny + 1
            global iDlny = iDlny + 1
            Dlny[iDlny, :] = [Layer[iLayer][i] Link[j]]
            nCLink[Layer[iLayer][i]] = nCLink[Layer[iLayer][i]]+1
            nCLink[Link[j]] = nCLink[Link[j]]+1
            # lines!(nuc[Dlny[iDlny,:],1], nuc[Dlny[iDlny,:],2])
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




Dlny = Dlny[1:iDlny,:]


# cell nuclei
scatter!(nuc[:,1], nuc[:,2], markersize = layerW/4.)
# Delaunay
for i in 1:iDlny
  lines!(nuc[Dlny[i,:],1], nuc[Dlny[i,:],2])
end
for i in 1:nLayers
     scatter!([nuc[Layer[i][1], 1]], [nuc[Layer[i][1],2]],
     markersize = .5, color = :red);
 end
display(s)
