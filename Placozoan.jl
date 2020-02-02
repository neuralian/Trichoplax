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
using Makie
using Colors

export Trichoplax,
                discworld, neighbour, makebody, findperimetervertices,
        smoothperimeter, makecellmap, makereceptivefields,
        drawdelaunaydisc, drawcells, drawskeleton, findperimeteredges,
        distalskeleton

struct Skeleton
  numlayers::Int64
  edgelength::Float64
  vertex::Array{Float64,2}
  edge::Array{Int64,2}
  layer::Vector{Array{Int64,1}}
  neighbour::Array{Int64}
end

struct Trichoplax
  numlayers::Int64                   # number of layers of cells
  celldiam::Float64
  skeleton::Skeleton
  vertex::Array{Float64}             # cell vertices
  cell::Array{Int64}                 # [i,j] index jth vertex of ith cell
  layer::Vector{Array{Int64,1} }     # [i][j] index jth cell in ith layer
  edge::Array{Int64}                 # [i,:] index links between cells
  skinvertex::Array{Int64}           #  skin vertex indices
  mapdepth::Int64                    # number of cell layers in sensory map
  mapvertex
  mapcell
end

normaldistribution = Normal()

numEdges(num_cells) = num_cells*6 + 9*num_cells*(num_cells-1)

function Skeleton(nlayers, edgelength)
    # Delaunay-triangulated disc
    # returns (vertex, edge, layer)
    # vertex = Float64 N x 2 coordinates
    # link =  Int64 Nx2  pairs of vertex indices
    # layer = ith row points to vertices in ith layer
    # MGP Dec 2019

    poswobble = .0025   # noise on cell location (re. edgewidth)
    nvertices = 3*nlayers*(nlayers-1)+1
    vertex = fill(0.0,nvertices , 2)
    edge = fill(0, 6*nvertices, 2)
    layer = fill(Int64[],nlayers)
    nlayer = Vector{Int64}(undef,nlayers)

    nlayer[1] = 1
    xs = 0.0
    for i in 2:nlayers
        fn = Float64(nlayer[i-1]) + 2.0π + xs
        n = round(fn)
        nlayer[i] = nlayer[i-1] + n
        xs = fn - n
    end

    nlayer = vcat([1], [6*i for i in 1:nlayers-1]) # n cells in layer

    nlinkpercell = fill(0, nvertices, 1)  # number of Delaunay links per cell
    icell = 1
    vertex[icell,:] = [0.0 0.0]  # cell at origin
    layer[1,1] = [1]
    iedge =0
    w = 0
    @inbounds for ilayer in 2:nlayers
        layer[ilayer,:] = [fill(0,nlayer[ilayer])]
        w = w + 2.0*(rand(1)[1] - 0.5)/nlayer[ilayer]
        @inbounds for j in 1:nlayer[ilayer]
            dw = rand(normaldistribution)
            while abs(dw)>3.0
                dw = rand(normaldistribution)
            end
            dw= .1*dw/(ilayer*edgelength)
            dr = rand(normaldistribution)
            while abs(dr)>3.0
                dr = rand(normaldistribution)
            end
            dr = poswobble*dr
            icell = icell + 1
            layer[ilayer][j] = icell
            vertex[icell,:] =
            (ilayer-1)*(1.0+dr)*edgelength.*
            [ cos(2π*((j-1)/nlayer[ilayer])+w+dw)
              sin(2π*((j-1)/nlayer[ilayer])+w+dw) ]

            # Delaunay triangulation
            if j>1
                iedge = iedge + 1
                edge[iedge, :] = [icell-1 icell]
                nlinkpercell[icell-1] = nlinkpercell[icell-1]+1
                nlinkpercell[icell] = nlinkpercell[icell]+1
            end
        end
        iedge = iedge + 1
        edge[iedge, :] = [icell layer[ilayer][1]]
        nlinkpercell[icell] = nlinkpercell[icell]+1
        nlinkpercell[layer[ilayer][1]] = nlinkpercell[layer[ilayer][1]]+1
    end

    # complete Delaunay triangles
    @inbounds for j in 1:6      # center cell special case, 6 links
        iedge = iedge + 1
        edge[iedge,:] = [1 j+1]
        nlinkpercell[1] = nlinkpercell[1]+1
        nlinkpercell[j+1] = nlinkpercell[j+1]+1
    end

    @inbounds for ilayer in 2:nlayers-1
        ilink = [ nlayer[ilayer+1] 1 2]
        link = layer[ilayer+1][ilink[1:3]]
        @inbounds for i in 1:length(layer[ilayer])
            #println((ilayer, link))
            @inbounds for j in 1:length(link)            #global iedge = iedge + 1
                iedge = iedge + 1
                edge[iedge, :] = [layer[ilayer][i] link[j]]
                nlinkpercell[layer[ilayer][i]] = nlinkpercell[layer[ilayer][i]]+1
                nlinkpercell[link[j]] = nlinkpercell[link[j]]+1
            end
            if i<length(layer[ilayer])
                link = link[end] .+ collect(0:(5-nlinkpercell[layer[ilayer][i+1]]))
            end
        end
    end
    edge = edge[1:iedge,:]

    neighbour=skeleton_neighbours(vertex,edge)

    #return (vertex, edge, layer)
    return Skeleton(nlayers, edgelength, vertex, edge, layer, neighbour)
end

function skeleton_neighbours(vertex, edge)
  # list neighbour vertices of each vertex
  # each vertex has 6 neighbours, except in outermost layer which have 3 or 4

  ncells = size(vertex, 1)
  neighbour = fill(0, ncells, 6)   # neighbour indices for each cell
  n_neighbour = fill(0, ncells)       # number of neighbour found

  # find neighbour of each cell
  @inbounds for i in 1:size(edge,1)
    a = edge[i,1]
    b = edge[i,2]
    n_neighbour[a] = n_neighbour[a]+1
    n_neighbour[b] = n_neighbour[b]+1
    neighbour[a, n_neighbour[a]] = b
    neighbour[b, n_neighbour[b]] = a
  end
  # sort into counterclockwise order
  # nb row 1 is already ccw
  @inbounds for i in 2:ncells
    n = count(x->(x>0), neighbour[i,:])
    dx = vertex[neighbour[i,1:n],1].-vertex[i,1]
    dy = vertex[neighbour[i,1:n],2].-vertex[i,2]
    theta = atan.(dy,dx)
    theta = theta .+ (π - atan(vertex[i,2], vertex[i,1]))
    i0 = findall(x->(x<-π), theta)
    theta[i0] = theta[i0] .+ 2π
    i1 = findall(x->(x>π), theta)
    theta[i1] = theta[i1] .- 2π
    order = sortperm(theta)
    neighbour[i,1:n] = neighbour[i, order]
  end
  neighbour
end

function distalskeleton(neighbour, layer)
    # skeleton triangles that define peripheral cell vertices

    # number of triangles is 2xnumber of vertices in next-to-last layer
    # plus the number of vertices with only 3 neighbour
    N = size(neighbour, 1)  # number of vertices in skeleton
    n3 = 0
    for i = 1:N
        if count(x -> (x > 0), neighbour[i, :]) == 3
            n3 = n3 + 1
        end
    end
    numtriangles = size(layer[end-1])[1] * 2 + n3
    triangle = fill(0, numtriangles, 3)

    # start at last skeleton vertex and work backwards
    n = N  # index skeleton vertex
    n_nbrs = count(x -> (x > 0), neighbour[n, :])
    i = 0   # index triangles
    while n_nbrs < 6     # outer ring of skeleton vertices have < 6 neighbour
        for j = 2:n_nbrs
            T = [n neighbour[n, j-1] neighbour[n, j]]
            alreadyfound = false
            for k = 1:i
                if any(isequal(T[1]), triangle[k, :]) &&
                   any(isequal(T[2]), triangle[k, :]) &&
                   any(isequal(T[3]), triangle[k, :])
                    alreadyfound = true
                    break
                end
            end
            if !alreadyfound
                i = i + 1
                triangle[i, :] = T
            end
        end
        n = n - 1
        n_nbrs = count(x -> (x > 0), neighbour[n, :])
    end
    triangle
end # distalskeleton

function makebody(skeleton)
  # construct cells
  # returns (vertex, cell, edge, layer)
  #  vertex = vertices of cell membrane (6 per cell)
  #  cell[i,:] = indices of 6 vertices of ith cell
  #  edge = Nx2 array of indices to ends of each membrane segment
  #  layer = index of cells in each layer
  #  e.g. cell[layer[end][1], :] points to vertices of first cell in last layer
  #  vertex[cell[layer[end][1], :],:]  is 6x2 array of vertices of this cell
  # MGP Dec 2019

  # tolerance for merging vertices.
  tol = 0.1*skeleton.edgelength

  ncells = skeleton.layer[end-1][end]       #
  vertex = fill(0.0, 6*ncells, 2 )  # vertex coordinates
  cell = fill(0, ncells, 6)         # vertex indices for each cell
  edge = fill(0, numEdges(ncells), 2)       # vertex indices for each edge

  nvertex = 0
  nedge = 0
  @inbounds for i in 1:ncells
      for j in 1:6
        k = skeleton.neighbour[i,j]
        m = skeleton.neighbour[i, j%6+1]
        x = sum(skeleton.vertex[[i k m],1])/3.
        y = sum(skeleton.vertex[[i k m],2])/3.
        vertex_exists = false
        @inbounds for n in 1:nvertex
            if sqrt( (x-vertex[n,1])^2 + (y-vertex[n,2])^2) < tol
                cell[i,j]=n  # point to existing vertex
                vertex_exists = true
                break
            end
        end
        if !vertex_exists # add new vertex and point to it
                nvertex = nvertex+1
                vertex[nvertex,:] = [x y]
                cell[i,j] = nvertex
        end # !vertex_exists

        # build list of edges (cell membrane segments)
        edge_exists = false
        if j>1
            e = sort([cell[i,j-1] cell[i,j]], dims=2)  # edge
            @inbounds for i1 in 1:nedge
                if vec(e)==edge[i1,:]
                    edge_exists = true
                    break
                end
            end # for i1
            if !edge_exists
            nedge = nedge + 1
            edge[nedge,:] = e
            end # !edge_exists
        end # if j>1
        # 6th edge links last and first vertex
        edge_exists = false
        if j==6
            e = sort([cell[i,j] cell[i,1]], dims=2)
            for i1 in 1:nedge
                if vec(e)==edge[i1,:]
                    edge_exists = true
                    break
                end
            end # for i1
            if !edge_exists
                nedge = nedge + 1
                edge[nedge,:] = e
            end # !edge_exists
        end # j==6
    end # for j
    end # for i
  (vertex[1:nvertex,:], cell, edge, skeleton.layer[1:(end-1)])
end #cells()

function findperimetervertices(cell, vlink, clayer)
    # find vertices on edge of disc
    n = length(clayer[end])
    outerperimeter = fill(0, n+6) # linked only to perimeter vertices
    cornerperimeter = fill(0, n)  # linked to 1 outer and 1 internal vertex
    innerperimeter = fill(0, n)   # internal vertices linked to corner vertices


    # outerperimeter (outer perimeter) vertices have no vlinks to internal vertices
    no = 0
    @inbounds for i in 1:n  # for each cell in outer layer
        @inbounds for v in 1:6  # vertices with < 3  links ..
            if sum(vlink[:]'.==cell[clayer[end][i],v])<3
                no = no + 1
                outerperimeter[no] = cell[clayer[end][i],v]
            end
        end
    end

    # cornerperimeter (inner perimeter) vertices
    # are vlinked to outerperimeter vertices + 1 internal vertex
    nc = 0
    @inbounds for i in 1:length(outerperimeter)
        @inbounds for j in 1:size(vlink,1) # for each vlink
            if  (vlink[j,1]==outerperimeter[i])    &  # 1st of jth vlink links to ith vertex
                 !any(outerperimeter.==vlink[j,2]) &   # and not to any other outerperimeter vertex
                 !any(cornerperimeter.==vlink[j,2])     # and hasn't already been found
               nc = nc+1
               cornerperimeter[nc] = vlink[j,2]         # found a corner perimeter vertex
           elseif (vlink[j,2]==outerperimeter[i]) &   # ditto for 2nd element in jth vlink
                 !any(outerperimeter.==vlink[j,1]) &
                 !any(cornerperimeter.==vlink[j,1])
               nc = nc+1
               cornerperimeter[nc] = vlink[j,1]
            end
        end
    end

    # internal vertices are vlinked to inner perimeter vertices
    # (used to align lateral membranes of perimeter cells orthogonal to skin )
    ni = 0
    @inbounds for i in 1:length(cornerperimeter)
        @inbounds for j in 1:size(vlink,1)
            if (vlink[j,1]==cornerperimeter[i])    & # 1st in jth vlink is a corner
                !any(outerperimeter.==vlink[j,2]) & # AND 2nd is not a perimeter cell
                !any(innerperimeter.==vlink[j,2])   # AND not already found
               ni = ni + 1
               innerperimeter[ni] = vlink[j,2]    # found interior link to corner
            elseif (vlink[j,2]==cornerperimeter[i])    & # 2nd in jth vlink is a corner
                   !any(outerperimeter.==vlink[j,1]) & # AND 1st is not a perimeter cell
                   !any(innerperimeter.==vlink[j,1])   # AND not already found
                  ni = ni + 1
                  innerperimeter[ni] = vlink[j,1]    # found interior link to corner
            end
        end
    end
    (outerperimeter, cornerperimeter, innerperimeter)
end

function smoothperimeter(vertex, op, cp, ip)
    # move perimeter vertices onto a circle

    No = length(op)
    Nc = length(ip)
    ro = fill(0.0, No)  # radii of outer vertices
    ri = fill(0.0, Nc)  # radii of internal vertices

    Ro = 0.0
    @inbounds for i in 1:No
        Ro  = Ro + sqrt(vertex[op[i],1]^2 + vertex[op[i],2]^2)
        ro[i] = sqrt(vertex[op[i],1]^2 + vertex[op[i],2]^2)
    end
    Rc = 0.0
    @inbounds for i in 1:Nc
        Rc = Rc + sqrt(vertex[cp[i],1]^2 + vertex[cp[i],2]^2)
        ri[i] = sqrt(vertex[ip[i],1]^2 + vertex[ip[i],2]^2)
    end

    # mean radius of perimeter vertices
    R = (Ro/No + Rc/Nc)/2.0

    # shift corner vertices to projection of inner vertices onto mean radius
    @inbounds for i in 1:Nc
        vertex[cp[i],1] = vertex[ip[i],1]*R/ri[i]
        vertex[cp[i],2] = vertex[ip[i],2]*R/ri[i]
    end

    @inbounds for i in 1:No
        vertex[op[i],1] = vertex[op[i],1]*R/ro[i]
        vertex[op[i],2] = vertex[op[i],2]*R/ro[i]
    end

    vertex
end

function Trichoplax(numlayers, celldiameter, mapdepth)

  skeleton = Skeleton(numlayers+1, celldiameter)

  # make celllayers layers of hexagonal cells (Voronoi tesselation)
  (vertex, cell, link, layer) = makebody(skeleton)

  # perimeter vertices in 3 classes
  perimetervertex = findperimetervertices(cell, link, layer)
  # ordered vertices on perimeter
  skinvertex = sort(vcat(perimetervertex[1],perimetervertex[2] ))

  (mapvertex,mapcell) = makecellmap(vertex, cell, layer, mapdepth);

  return Trichoplax(numlayers, celldiameter, skeleton,
       vertex, cell, layer, link,
       skinvertex, mapdepth, mapvertex, mapcell)
end

function makecellmap(vertex, cell, clayer, mapwidth)
    # reflect outer layers of cells into environment
    # by reflecting mapWidth layers of cells around body perimeter
    # into annular region of width senseRange

    ncells = size(cell,1)

    #number of cell layers in body
    nlayers = Int64((3+sqrt(9+12*(ncells-1)))/6)

    #number of cells in map (outer celldepth layers)
    nmap = ncells - (3*(nlayers-mapwidth)*(nlayers-mapwidth-1) + 1)

    # duplicate body
    mapcell = fill(0, nmap, 6)

    nvertex = size(vertex,1)
    @inbounds for i in 1:nmap  #size(cell,1)
        @inbounds for j in 1:6
            mapcell[i,j] = cell[ncells-i+1,j]
        end
    end

    # minimum cell index in pCell
    i0 = findmin(mapcell)[1]

    # copy the required vertices, shift indices to coincide
    #pVertex = copy(vertex[i0:end,:])
    nmapvertex = nvertex-i0+1
    mapvertex = fill(0.0, nmapvertex,2 )
    mapcell = mapcell.-(i0-1)

    # project vertices beyond skin
    Rs = sqrt(vertex[end,1]^2 + vertex[end,2]^2) # skin radius
    Rc = sqrt(vertex[2,1]^2 + vertex[2,2]^2)  # cell radius
    @inbounds for i in 1:nmapvertex
       R0 = sqrt(vertex[i0+i-1,1]^2 + vertex[i0+i-1,2]^2)
       R1 = Rs + (Rs-R0)*(Rs/R0)^2
       mapvertex[i,:] = vertex[i0+i-1,:].*R1/R0
    end

    (mapvertex, mapcell)
end

function makereceptivefields(vertex, cell, clayer, mapwidth)
    # construct receptive field centers
    # by reflecting centroids of mapping cells
    # into annular region of width senseRange

    ncells = size(cell,1)

    #number of cell layers in body
    nlayers = Int64((3+sqrt(9+12*(ncells-1)))/6)

    # number of cells in map (outer celldepth layers)
    nmap = ncells - (3*(nlayers-mapwidth)*(nlayers-mapwidth-1) + 1)

    # index of cell (centroid) preceding first cell in map
    i0 = 3*(nlayers-mapwidth)*(nlayers-mapwidth-1)+1

    # array of rf centers
    rfcenter = fill(0.0, nmap,2)

    # project vertices beyond skin
    Rs = sqrt(vertex[end,1]^2 + vertex[end,2]^2) # skin radius
    Rc = sqrt(vertex[2,1]^2 + vertex[2,2]^2)  # cell radius
    @inbounds for i in 1:nmap
        x = mean(vertex[cell[i0+i,:],1])
        y = mean(vertex[cell[i0+i,:],2])
        R0 = sqrt(x^2 + y^2)
        R1 = Rs + (Rs-R0)*(Rs/R0)^2
       rfcenter[i,:] = [x*R1/R0 y*R1/R0]
    end

    return rfcenter
end

# function drawdelaunaydisc(dlink,color = RGB(0.25, 0.25, 0.25), linewidth = 0.5)
#     @inbounds for i = 1:size(dlink, 1)
#       lines!(
#         cellcentroid[dlink[i, :], 1],
#         cellcentroid[dlink[i, :], 2],
#         color = color, linewidth = linewidth
#         )
#     end
# end



# function drawcells(vertex, vlink)
#     # draw links between cell vertices (cytoskeleton)
#     for i in 1:size(vlink, 1)
#         lines!(vertex[vlink[i,:],1], vertex[vlink[i,:],2])
#     end
# end

end  # module Placozoan
