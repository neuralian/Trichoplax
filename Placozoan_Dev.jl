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

function cells(nucleus, neighbour, dlayer)
  # construct cells
  # returns (vertex, cell, vlink, clayer)
  #  vertex = vertices of cell membrane (6 per cell)
  #  cell[i,:] = indices of 6 vertices of ith cell
  #  vlink = Nx2 array of indices to ends of each membrane segment
  #  clayer = index of cells in each layer
  #  e.g. cell[clayer[end][1], :] points to vertices of first cell in last layer
  #  vertex[cell[clayer[end][1], :],:]  is 6x2 array of vertices of this cell
  # MGP Dec 2019

  # tolerance for merging vertices.
  # nb assuming nucleus_1 is next to nucleus_2, so distance between them
  # is cell diameter
  tol = 0.1*sqrt((nucleus[1,1]-nucleus[2,1])^2 + (nucleus[1,2]-nucleus[2,2])^2)

  nCells = size(nucleus,1) - length(dlayer[end])
  vertex = fill(0.0, 6*nCells, 2 )  # vertex coordinates
  cell = fill(0, nCells, 6)         # vertex indices for each cell
  vlink = fill(0, 6*nCells, 2)       # vertex indices for each vlink

  nVertex = 0
  nvlink = 0
  for i in 1:nCells
      for j in 1:6
        k = neighbour[i,j]
        m = neighbour[i, j%6+1]
        x = (nucleus[i,1] + nucleus[k,1] + nucleus[m,1])/3.
        y = (nucleus[i,2] + nucleus[k,2] + nucleus[m,2])/3.
        vertexExists = false
        for n in 1:nVertex
            if sqrt( (x-vertex[n,1])^2 + (y-vertex[n,2])^2) < tol
                cell[i,j]=n  # point to existing vertex
                vertexExists = true
                break
            end
        end
        if !vertexExists # add new vertex and point to it
                nVertex = nVertex+1
                vertex[nVertex,:] = [x y]
                cell[i,j] = nVertex
        end # !vertexExists

        # build list of vlinks (cell membrane segments)
        vlinkExists = false
        if j>1
            e = sort([cell[i,j-1] cell[i,j]], dims=2)  # vlink
            for i1 in 1:nvlink
                if vec(e)==vlink[i1,:]
                    vlinkExists = true
                    break
                end
            end # for i1
            if !vlinkExists
            nvlink = nvlink + 1
            vlink[nvlink,:] = e
            end # !vlinkExists
        end # if j>1
        # 6th vlink links last and first vertex
        vlinkExists = false
        if j==6
            e = sort([cell[i,j] cell[i,1]], dims=2)
            for i1 in 1:nvlink
                if vec(e)==vlink[i1,:]
                    vlinkExists = true
                    break
                end
            end # for i1
            if !vlinkExists
                nvlink = nvlink + 1
                vlink[nvlink,:] = e
            end # !vlinkExists
        end # j==6
    end # for j
    end # for i
  (vertex[1:nVertex,:], cell, vlink[1:nvlink, :], dlayer[1:(end-1)])
end #cells()


function find_perimeter(cell, vlink, clayer)
    # find vertices and vlinks on edge of disc
    operimeter = fill(0, 2*size(clayer[end],1)) # outer
    iperimeter = fill(0, 2*size(clayer[end],1))# inner (linked to internal vertex)
    no = 0
    for i in 1:length(clayer[end])  # for each cell in outer layer
        for v in 1:6  # find vertices with fewer than 3 incoming links ..
            if sum(vlink[:]'.==cell[clayer[end][i],v])<3
                no = no + 1
                operimeter[no] = cell[clayer[end][i],v]
            end
        end
    end
    # iperimeter vertices are vlinked to operimeter vertices + 1 internal vertex
    ni = 0
    for i in 1:length(operimeter)
        for j in 1:size(vlink,1)
            if (vlink[j,1]==operimeter[i]) & !any(operimeter.==vlink[j,2])
                ni = ni+1
                iperimeter[ni] = vlink[j,2]
            elseif (vlink[j,2]==operimeter[i]) & !any(operimeter.==vlink[j,1])
                ni = ni+1
                iperimeter[ni] = vlink[j,1]
            end
        end
    end

    (operimeter[1:no], unique(iperimeter[1:ni]))
end

function round_perimeter(vertex, op, ip)
    # move perimeter vertices onto a circle

    No = length(op)
    Ni = length(ip)
    ri = fill(0.0, Ni)
    ro = fill(0.0, No)
    for i in 1:No
        ro[i] = sqrt(vertex[op[i],1]^2 + vertex[op[i],2]^2)
    end
    for i in 1:Ni
        ri[i] = sqrt(vertex[ip[i],1]^2 + vertex[ip[i],2]^2)
    end

    # average radius of perimeter vertices
    r =  (sum(ro)+sum(ri))/(No+Ni)

    # shift inner vertices to mean radius
    for i in 1:Ni
        vertex[ip[i],1] = vertex[ip[i],1]*r/ri[i]
        vertex[ip[i],2] = vertex[ip[i],2]*r/ri[i]
    end

    for i in 1:No
        vertex[op[i],1] = vertex[op[i],1]*r/ro[i]
        vertex[op[i],2] = vertex[op[i],2]*r/ro[i]
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

nLayer = 25
layerWidth = 5.0
print("DelaunayDisc:")
@time DD = delaunayDisc(nLayer, layerWidth)
print("neighbours:")
@time nbrs=neighbours(DD...)


cellNucleus = DD[1];
dlink = DD[2];
dlayer = DD[3];
print("cells:")
@time C = cells(cellNucleus, nbrs, dlayer)
vertex = C[1]
cell = C[2]
vlink = C[3]
clayer = C[4]

# Draw
s = Scene(resolution = (800,800), scale_plot = false)

# draw cell nuclei
scatter!(
  cellNucleus[1:size(cell,1), 1],
  cellNucleus[1:size(cell,1), 2],
  markersize = layerWidth / 6.0,
  color = RGB(1.0, .5, 0.2),
)

function draw_cells(vertex, cell)
    @inbounds for i in 1:size(cell,1)
        lines!(vertex[cell[i,[1:6; 1]],1], vertex[cell[i,[1:6; 1]],2],
                color = RGB(.65, .45, 0.0))
    end
end



function draw_vlinks(vertex, vlink)
    for i in 1:size(vlink, 1)
        lines!(vertex[vlink[i,:],1], vertex[vlink[i,:],2])
    end
end

p = find_perimeter(cell, vlink, clayer)
op = p[1]
ip = p[2]

vertex = round_perimeter(vertex, op, ip)


draw_cells(vertex,cell)

# scatter!(vertex[op,1], vertex[op,2], markersize = layerWidth/4., color = :red)
# scatter!(vertex[ip,1], vertex[ip,2], markersize = layerWidth/4., color = :green)
display(s)
