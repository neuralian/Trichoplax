using Makie


nLayers = 3
layerW = 5.0
nCells = Int64(round(π*(nLayers+1)^2))
nuc = fill(0.0,nCells , 2)
Dlny = fill(0, 6*nCells, 2)
Layer = fill(Int64[],nLayers+1)

iCell = 1
nuc[iCell,:] = [0.0 0.0]  # cell at origin
Layer[1,1] = [1]
iDlny =0
for iLayer in 1:nLayers
    nLayer = Int64(round(2π*iLayer))
    Layer[iLayer+1,:] = [fill(0,nLayer)]
    w = rand(1)[1]/2
    for j in 1:nLayer
        global iCell = iCell + 1
        Layer[iLayer+1][j] = iCell
        nuc[iCell,:] =
        iLayer*layerW.*[cos(2π*((j+w)/nLayer))
                       sin(2π*((j+w)/nLayer))]

        # Delaunay triangulation
        if j>1
            global iDlny = iDlny + 1
            Dlny[iDlny, :] = [iCell-1 iCell]
        end

    end
end
nuc = nuc[1:iCell,:]
Dlny = Dlny[1:iDlny,:]

s = Scene(scale_plot = false)
# cell nuclei
scatter!(nuc[:,1], nuc[:,2], markersize = layerW/4.)
# delaunay
for i in 1:iDlny
  lines!(nuc[Dlny[i,:],1], nuc[Dlny[i,:],2])
end

display(s)
