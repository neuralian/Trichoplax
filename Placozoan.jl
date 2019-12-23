# Placozoan.jl Package/module
#  To make the package available:
#   At REPL prompt: julia> push!(LOAD_PATH, <path_to_package_folder>)
# In Juno:
#    Right-click on folder containing the package.
#    Select "Juno: Work in folder"
#    At REPL prompt:  julia> push!(LOAD_PATH, pwd())
#
# MGP 2019-2020

module Placozoan
using Distributions

export delaunayDisc

Ndist = Normal()

function delaunayDisc(nLayers, layerWidth)
    # returns (nucleus, dlink, layer)
    # nucleus = nCells x 2 coordinates of ith cell nucleus
    # dlink = ndlink x indices of cell nuclei defining the dlink
    # layer = ith row indexes cell nuclei in each layer
    # MGP Dec 2019
a = 0.005    # noise on cell location (re. layerWidth)
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
    w = w + 2π*(rand(1)[1] - 0.5)/nLayer[iLayer]
    @inbounds for j in 1:nLayer[iLayer]
        dr = rand(Ndist)
        while abs(dr)>3.0
            dr = rand(Ndist)
        end
        dw = rand(Ndist)
        while abs(dw)>3.0
            dw = rand(Ndist)
        end
        w = w + a*dw
        iCell = iCell + 1
        Layer[iLayer][j] = iCell
        nucleus[iCell,:] =
        (iLayer-1)*(1.0+a*dr)*layerWidth.*
        [cos(2π*((j-1)/nLayer[iLayer])+w)
                       sin(2π*((j-1)/nLayer[iLayer])+w)]

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


end  # module Placozoan
