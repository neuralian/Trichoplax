# Placozoan.jl Package/module
#  To make the package available:
#   At REPL prompt: julia> push!(LOAD_PATH, <path_to_package_folder>)
# In Juno:
#    Right-click on folder containing the package.
#    Select "Juno: Work in folder"
#    At REPL prompt:  julia> push!(LOAD_PATH, pwd())
#
# MGP 2019-2020

# module Placozoan
#
# using Distributions
using Makie
using Colors

# export Trichoplax,
#                 discworld, neighbour, makebody, findperimetervertices,
#         smoothperimeter, makecellmap, makereceptivefields,
#         drawdelaunaydisc, drawcells, drawskeleton, findperimeteredges,
#         distalskeleton


struct Skeleton
    # Delaunay-triangulated disc
    # triangles define vertices of hexagonal cells (Voronoi tesselation)
    vertex::Array{Float64,2}          # ?x2 vertices
    link::Array{Int64,2}              # ?x2 edges forming Delaunay Δ-ation
    edgelength::Array{Float64,1}      # nominal edge length (= cell diameter)
    layercount::Array{Int64,1}        # vertices in ith layer
    neighbour::Array{Int64,2}         # 6 neighbours of each internal vertex
end

struct Trichoplax
    nlayers::Int64
    celldiameter::Float64      # nominal cell diameter when constructed
    # skeleton::Skeleton
    triangle::Array{Int64,2}   # skeleton Delaunay triangles
    skintriangle::Array{Int64,2} # triangles for skin vertices
    vertex::Array{Float64,2}   # cell vertices, each is centroid of triangle
    cell::Array{Int64, 2}      # 6 vertis for each cell
    edge::Array{Int64}         # [i,:] index links between cells
    edgelength::Array{Float64,1}  # edge rest lengths
    volume::Array{Float64,1}   # volume (area) of each cell
    # skinvertex::Array{Int64}           #  skin vertex indices
    # mapdepth::Int64                    # number of cell layers in sensory map
    # mapvertex
    # mapcell
    skin::Array{Int64}         # index to cell vertices on exterior surface
    n_edges2::Array{Int64,1}   # number of edges at ith vertex
    edge2::Array{Int64,2}     # index of edges at ith vertex
    n_neighbour::Array{Int64,1}         # number of neighbours for each vertex
    neighbour::Array{Int64,2}  # neighbours of each vertex
    n_vertexcells::Array{Int64,1}  # number of cells containing each vertex
    vertexcells::Array{Int64,2} # index of cells containing each vertex
    skin_neighbour::Array{Int64,2}  # 2 skin neighbours for each skin vertex
    k2::Array{Float64,1}  # half of cytoskeleton spring constant (k/2)
    ρ::Array{Float64,1}   # cell pressure constant \rho (energy/volume)
    σ::Array{Float64,1}   # surface energy density
    # nb fields are immutable
    #    Declaring k2 & σ as arrays allows access via e.g. trichoplax.k2[]
end

# normaldistribution = Normal()

#numLinks(num_cells) = num_cells*6 + 9*num_cells*(num_cells-1)

function Skeleton(n_cell_layers, cell_diameter)

    layercount = skeletonlayercount(n_cell_layers)
    vertex = skeletonvertices(layercount, cell_diameter)
    neighbour = skeletonvertexneighbours(vertex, layercount)
    link = links(vertex,neighbour, layercount)
    # distalΔ = distalskeleton(nbr, layer)

    return Skeleton(vertex, link, [cell_diameter], layercount, neighbour)
end

function Trichoplax(n_cell_layers, cell_diameter)

    # skeleton is a triangulated disc  (Delaunay triangulation)
    # triangle edge length = cell_diameter
    # each internal vertex (ie except for outer layer) will be a cell centre
    skeleton = Skeleton(n_cell_layers, cell_diameter)

    #  vertices of hexagonal cells are at centres of triangles
    # (Voronoi tesselation)
    # cell[i,:] indexes the vertices of ith cell
    # triangle[i,:] indexes skeleton vertices for rapid updating cell vertices
    #               when skeleton moves
   (vertex, cell, triangle, skintriangle) = makecells(skeleton)

  edge = celledges(cell)
  edgelength = edgelengths(edge, vertex)

  nVertices = size(vertex,1)
  (n_neighbour, neighbour) = cellvertexneighbours(nVertices, edge)
  (n_vertexcells, vertexcells) = cellscontainingcellvertices(nVertices, cell)

  (n_edges2, edge2) = edges2vertex(nVertices, edge)

  (skin, skinvertex) =  getskin(vertex, skeleton.vertex, skintriangle)

  skin_neighbour = skinneighboursofskinvertices(neighbour, skin)

  volume = ones(size(cell,1))*cellvolume(vertex, cell)[1]

  skeleton_springconstant = 2.0
  cell_surface_energy_density = 25.0
  cell_pressureconstant = 1.0

  return Trichoplax(    n_cell_layers,
                        cell_diameter,
                        triangle,
                        skintriangle,
                        vertex,
                        cell,
                        edge,
                        edgelength,
                        volume,
                        skin,
                        n_edges2,
                        edge2,
                        n_neighbour,
                        neighbour,
                        n_vertexcells,
                        vertexcells,
                        skin_neighbour,
                        [skeleton_springconstant/2.0],
                        [cell_pressureconstant],
                        [cell_surface_energy_density]

                        )
end

function skeletonlayercount(n_cell_layers)
    # number of vertices in each skeleton layer given number of cell layers

   return vcat([1], [6*i for i in 1:n_cell_layers])
end

function skeletonvertices(layercount, edgelength)
    #
    nlayers = length(layercount)
    nvertices = sum(layercount)
    vertex = fill(0.0, nvertices, 2)  # vertex 1 is [0,0]

    i = 1                         # initial vertex index
    for j in 2:nlayers            # for each layer
        r = (sqrt(3)/2)*(j-1)*edgelength      # base radius
        for k in 1:layercount[j]  # for each vertex in layer
            i = i + 1
            θ = 2π*(k-1)/layercount[j]
            vertex[i,:] = r*[cos(θ) sin(θ)]
        end
    end
    return vertex
end

function vertexdistance(v)
    # D[i,j] = distance from ith to jth vertex

    n = size(v,1)
    D = fill(0.0, n, n)
    for i in 1:n
        for j in 1:n
            D[i,j] = sqrt( (v[i,1]-v[j,1])^2 + (v[i,2]-v[j,2])^2 )
        end
    end

    return D
end

function skeletonvertexneighbours(v, layercount)
    # 6 neighbours of each skeleton vertex, not including outer layer
    # sorted into counterclockwise order

    n = sum(layercount[1:(end-1)])   # number of vertices with 6 neighbours
    D = vertexdistance(v)
    neighbour = fill(0, n, 6)

    # closest 6 neighbours of each vertex
    for i in 1:n
        order = sortperm(D[:,i])  # list vertices in order of distance from ith
        neighbour[i, :] = order[2:7] # nb closest vertex order[1]==self
    end

    # sort neighbours anticlockwise around self
    for i in 1:n
      dx = v[neighbour[i,:],1].-v[i,1]
      dy = v[neighbour[i,:],2].-v[i,2]
      θ = atan.(dy,dx)
      θ = θ .+ (π - atan(v[i,2], v[i,1]))
      i0 = findall(x->(x<-π), θ)
      θ[i0] = θ[i0] .+ 2π
      i1 = findall(x->(x>π), θ)
      θ[i1] = θ[i1] .- 2π
      order = sortperm(θ)
      neighbour[i,:] = neighbour[i, order]
    end

    return neighbour
end

function links(v,neighbour, layercount)
    # links between skeleton vertices (= links between cells)

    numLinks(nlayers) = 3*nlayers*(3*nlayers-1)+layercount[end]
    nCells = sum(layercount[1:(end-1)])
    nLinks = numLinks(size(layercount,1)-1)

    link = fill(0, nLinks, 2)
    countLink = 0  # number of links found

    for i in 1:nCells
        for j in 1:6
            # candidate link from ith cell to its jth neighbour
            candidate = sort([i neighbour[i,j]], dims=2)
            already_found = false
            for k in 1:countLink
                if candidate==link[k,:]'
                    already_found = true
                    break
                end
            end
            if !already_found
                countLink = countLink + 1
                link[countLink,:] = candidate
            end
        end
    end

    # add peripheral ring
    i1 = sum(layercount)
    i0 = i1 - layercount[end]+1
    countLink = countLink + 1
    link[countLink, :] = [i1 i0]
    for i in (i0+1):i1
        countLink = countLink + 1
        link[countLink,:] = [i i-1]
    end


   return link
end

function cellvertexneighbours(nVertex, edge)
    # list neighbours of each vertex (vertices with an edge to this)

    #nVertex = size(vertex, 1)
    nEdge = size(edge,1)
    println(nVertex, ", ", nEdge)
    neighbour = fill(0, nVertex, 3)
    n_neighbour = fill(0, nVertex)  # number of edges connected to ith vertex
    for i in 1:nVertex
        neighbourcount = 0
        for j in 1:nEdge
            if any(edge[j,:].==i)      # found edge containing ith vertex
                # next neighbour is the other vertex in the link
                # (the one that is not i)
                neighbourcount = neighbourcount + 1
                neighbour[i, neighbourcount] = edge[j,findfirst(edge[j,:].!=i)]
            end
        end
        n_neighbour[i] = neighbourcount
    end
    return (n_neighbour, neighbour)
end


function edges2vertex(nVertex, edge)
    # list edges that connect to each vertex

    nEdge = size(edge,1)
    edge2 = fill(0, nVertex, 3)
    n_edges2 = fill(0, nVertex)  # number of edges connected to ith vertex
    for i in 1:nVertex
        edgecount = 0
        for j in 1:nEdge
            if any(edge[j,:].==i)      # jth edge contains ith vertex
                edgecount = edgecount + 1
                edge2[i, edgecount] = j
            end
        end
        n_edges2[i] = edgecount
    end
    return (n_edges2, edge2)
end

function cellscontainingcellvertices(nVertex, cell)
    # list cells containing each vertex

    #nVertex = size(vertex, 1)
    nCell = size(cell,1)
    cellshere = fill(0, nVertex, 3)
    n_cellshere = fill(0, nVertex)
    for i in 1:nVertex
        cellcount = 0
        for j in 1:nCell
            if any(cell[j,:].==i)
                cellcount = cellcount + 1
                cellshere[i, cellcount] = j
            end
        end
        n_cellshere[i] = cellcount
    end
    return (n_cellshere, cellshere)
end

function skinneighboursofskinvertices(neighbour, skin)
    # index skin vertices connected to each skin vertex
    # neighbour = index of neighbours of each vertex
    #   ( output from cellvertexneighbours() )
    # skin = index of vertices in skin
    #   ( output from getskin() )

    # nb inefficient code but it's simpler to search everywhere than only
    #    where the answer could be

    nSkinvertices = size(skin,1)
    skin_neighbour = fill(0, nSkinvertices, 2) # each skin vertex has 2 skin neighbours
    for i in 1:nSkinvertices
        n_neighbours = 0
        for j in 1:3
            if any(skin.==neighbour[skin[i],j])
                n_neighbours = n_neighbours + 1
                skin_neighbour[i,n_neighbours] = neighbour[skin[i],j]
            end
        end
    end
    return skin_neighbour
end

function cellvolume(vertex, cell)
    # cell volumes (area in 2D) from cell vertices

    nCells = size(cell,1)
    volume = fill(0.0, nCells )
    for i in 1:nCells
        v = 0.0
        for j in 1:6
            k = j % 6 + 1
            v = v + vertex[cell[i,j],1]*vertex[cell[i,k],2] -
                    vertex[cell[i,j],2]*vertex[cell[i,k],1]
        end
        volume[i] = abs(v/2.0)
    end

    return volume
end

function cellvolume(vertex)
    # cell volume given 6 vertices

    v = 0.0
    for j in 1:6
        k = j % 6 + 1
        v = v + vertex[j,1]*vertex[k,2] - vertex[k,1]*vertex[j,2]
    end

    return abs(v/2.0)
end

function makecells(skeleton)
    # cell vertices given skeleton vertices v & neighbours nbr
    # also returns index of vertices for each cell
    #  and index of skeleton triangles containing cell vertices
    #  so that cell vertices can be quickly updated when skeleton moves

    tol = .01*skeleton.edgelength[]

    nlayers = size(skeleton.layercount,1)-1

    nCells = size(skeleton.neighbour,1)
    nVertices = 6*nlayers^2
    vertex = fill(0.0, nVertices, 2 )
    triangle = fill(0, nVertices, 3) # skeleton triangle containing vertex
    cell = fill(0, nCells, 6)  # index 6 vertices per cell

    v = skeleton.vertex  # alias for code clarity (nb. no copy, v is a pointer)
    nbr = skeleton.neighbour

    n_vertex = 0
    n_cell = 0
    found_index = 0
    for i_cell in 1:nCells
        for j in 1:6
            k = j==6 ? 1 : j+1
            # compute cell vertex coords at centre of skeleton (i,j,k) triangle
            x = sum([v[i_cell,1] v[nbr[i_cell,j],1] v[nbr[i_cell,k],1]])/3.
            y = sum([v[i_cell,2] v[nbr[i_cell,j],2] v[nbr[i_cell,k],2]])/3.

            # is vertex already found?
            already_found = false
            for i_vertex in 1:n_vertex
                if (abs(x-vertex[i_vertex,1])<tol) &&
                   (abs(y-vertex[i_vertex,2])<tol)
                    already_found = true
                    found_index = i_vertex
                    break
                end
            end
            if already_found                   # point to existing vertex
                cell[i_cell,j] = found_index
            else                               # create new vertex
                n_vertex = n_vertex + 1        # and record triangle containing it
                vertex[n_vertex,:] = [x,y]
                cell[i_cell, j] = n_vertex
                triangle[n_vertex,:] = [i_cell nbr[i_cell,j] nbr[i_cell,k]]
            end

        end
    end

     # sort triangles that define skin vertices into anticlockwise order
     skinstart = 6*(nlayers-1)^2+1 # index to first skin triangle

     # extract outer layer of triangles (which define skin vertices)
     # in anticlockwise order
     ii = skinstart:nVertices
     θ = atan.(vertex[ii, 2], vertex[ii,1])
     skintriangle = triangle[ii[sortperm(θ)],:]

    return (vertex, cell, triangle, skintriangle)

end

function celledges(cell)
    # construct cytoskeleton links (cell edges) from cell vertex indices
    # cell  = nCell x 6 index of vertices for each cell

    nCells = size(cell,1)
    nEdges = nCells*6      # upper bound for now
    edge = fill(0, nEdges, 2)
    nfoundedges = 0

    for i in 1:nCells
        e0 = cell[i,end]
        for j in 1:6
            e1 = cell[i,j]
            candidate_edge = sort([e0 e1], dims=2)
            already_found = false
            for k in 1:nfoundedges
                if candidate_edge == edge[k,:]'
                    already_found = true
                    break
                end
            end
            if !already_found
                nfoundedges = nfoundedges + 1
                edge[nfoundedges,:] = candidate_edge
            end
            e0 = e1
        end
    end

    return edge[1:nfoundedges,:]
end

function edgelengths(edge::Array{Int64,2}, vertex::Array{Float64,2})

    nedge = size(edge,1)
    edgelength = fill(0.0, nedge)
    for i in 1:nedge
        edgelength[i] = sqrt( (vertex[edge[i,1],1] - vertex[edge[i,2],1])^2 +
                              (vertex[edge[i,1],2] - vertex[edge[i,2],2])^2 )
    end
    return edgelength
end

function cellverticesfromskeleton(trichoplax)
    #    (re-)compute cell vertices from skeleton vertices

    n = size(trichoplax.triangle, 1)
    vertex = fill(0.0, n, 2)
    v = trichoplax.skeleton.vertex
    for i in 1:n
        trichoplax.vertex[i,:] = sum(v[trichoplax.triangle[i,:],:], dims=1)/3.
    end

    return trichoplax
end


function getskin(cellvertex, skeletonvertex, skintriangle)
    # compute external vertex coordinates (skin vertices)
    # from skeleton skin triangles
    # and index these vertices in cell vertex array

    n = size(skintriangle, 1)
    skinvertex = fill(0.0, n, 2) # vertex coords
    skin = fill(0,n)             # index skin vertices in trichoplax.vertex
    nfoundskin = 0
    vc = cellvertex
    vs = skeletonvertex
    # vertex match tolerance 1% of distance between vertex 1 and 2
    tol = ( (vc[1,1] - vc[2,1])^2 + (vc[1,2] - vc[2,2])^2 )*(.01)^2
    for i in 1:n
        skinvertex[i,:] = sum(vs[skintriangle[i,:],:], dims=1)/3.
        for j in 1:size(vc,1)
            if ((skinvertex[i,1] - vc[j,1])^2 +
                (skinvertex[i,2] - vc[j,2])^2 ) < tol
                nfoundskin = nfoundskin + 1
                skin[nfoundskin] = j
                break
            end
        end
    end

    return (skin, skinvertex)
end

function draw(trichoplax::Trichoplax, color=:black, linewidth = .25)

    @inbounds for i in 1:size(trichoplax.cell,1)
        lines!(trichoplax.vertex[trichoplax.cell[i,[1:6; 1]],1],
               trichoplax.vertex[trichoplax.cell[i,[1:6; 1]],2],
                color = color, linewidth=linewidth)
    end
end

function drawskeleton(skeleton::Skeleton,
         color = RGB(.25,.65,.25), linewidth = 0.25)

  for i in 1:size(link,1)
      lines!(   skeleton.vertex[skeleton.link[i, :],1],
                skeleton.vertex[skeleton.link[i, :],2],
                color=color, linewidth=linewidth)
  end
end

function plotneighbours(v, nbrs)
    for i in 1:size(nbrs,1)
        for j = 1:6
            plot!(vec([v[i,1] v[nbr[i,j],1]]), vec([v[i,2] v[nbr[i,j],2]]))
        end
    end
end

function skeletonEnergy(v::Array{Float64,2}, trichoplax::Trichoplax)
    #  potential energy as a function of trichoplax shape (vertex coords)
    #   = elastic energy in skeleton deformation + surface energy.
    #   skeleton edges are linear springs Es = (1/2)k(r-ro)^2
    #   surface energy σ per unit length of exposed membrane
    #
    #  vertex is skeleton vertex coords [x1 x2 ... x_nv y1 y2 ... y_nv]
    #   i.e. nv X 2 vertex array flattened to a vector,
    #        vec(trichoplax.skeleton.vertex)

    # elastic energy in skeleton
    Es = 0.0
    link = trichoplax.skeleton.link   # alias for code readability
    nVertex = size(trichoplax.skeleton.vertex,1)
    p = 8
    q = 1.0/p
    for i in 1:size(trichoplax.skeleton.link,1)
        r = (  ( v[link[i,1],1] - v[link[i,2],1] ) ^p +
                   ( v[link[i,1],2] - v[link[i,2],2] ) ^p )^q
        Es = Es + trichoplax.k2[]*(r - trichoplax.skeleton.edgelength[])^2
    end

    # surface energy of external membranes
    Δ = trichoplax.skintriangle
    Lx = 0.0
    v0 = sum(v[Δ[end,:], :], dims=1)/3.0  # last skin vertex
    for i in 1:size(Δ,1)
        v1 = sum(v[Δ[i,:], :], dims=1)/3.0  # next skin vertex
        Lx = Lx + sqrt( (v1[1]-v0[1])^2 + (v1[2]-v0[2])^2)
        v0 = v1
    end

    # println(Es, ", ", trichoplax.σ[]*Lx)
     return (Es +  trichoplax.σ[]*Lx)
end

function shapeEnergy(trichoplax::Trichoplax)
    #  potential energy as a function of trichoplax shape (cell vertex coords)
    #   = elastic energy in cytoskeleton + turgor pressure + surface energy.
    #   skeleton edges are linear springs Es = (1/2)k(r-ro)^2
    #   surface energy σ per unit length of exposed membrane

    # update cell vertices from skeleton vertices
    #trichoplax = cellverticesfromskeleton(trichoplax)

    # elastic energy in cytoskeleton
    # TODO: using skeleton edge length as proxy for cell edge length
    Es = 0.0
    v = trichoplax.vertex
    edge = trichoplax.edge
    @inbounds for i in 1:size(edge,1)
        r = sqrt(  ( v[edge[i,1],1] - v[edge[i,2],1] ) ^2 +
                   ( v[edge[i,1],2] - v[edge[i,2],2] ) ^2 )
        Es = Es + trichoplax.k2[]*(r - trichoplax.skeleton.edgelength[])^2
    end

    # pressure
    Ep = 0.0
    for i in 1:size(trichoplax.cell, 1)
        Ep = Ep + (cellvolume(v,i)-trichoplax.volume[i])^12
    end


    # surface energy of external membranes
    Lx = 0.0
    skin = trichoplax.skin
    v0 = v[skin[end],:]
    @inbounds for i in 1:length(skin)
        v1 = v[skin[i],:]
        Lx = Lx + sqrt( (v1[1]-v0[1])^2 + (v1[2]-v0[2])^2)
        v0 = v1
    end

    # println(Es, ", ", trichoplax.σ[]*Lx)
     return (Es +  trichoplax.ρ[]*Ep + trichoplax.σ[]*Lx)
end

function shapeEnergyGradient(dv::Float64, trichoplax::Trichoplax)

    v = trichoplax.vertex
    n = size(v)
    ∇E = fill(0.0, n )



    E0 = shapeEnergy(trichoplax)
    dv = dv*trichoplax.skeleton.edgelength[]
    @inbounds for i in 1:n[1]
        @inbounds for j in 1:n[2]
            v[i,j] = v[i,j] + dv  # + perturb (i,j)th coordinate
            ∂Eplus = shapeEnergy(trichoplax)
            v[i,j] = v[i,j] - 2.0*dv  # - perturb (i,j)th coordinate
            ∂Eminus = shapeEnergy(trichoplax)
            ∇E[i,j] = (∂Eplus-∂Eminus)/(2.0*dv) # ∂E/∂v_ij
            v[i,j] = v[i,j] + dv      # put (i,j)th coordinate back
        end
    end
    ∇E
end

function localShapeEnergy(i::Int64, trichoplax::Trichoplax)
    #  component of potential energy that depends on the ith vertex

    Es = 0.0   # elastic energy in edges
    Le = 0.0   # edge length
    v = trichoplax.vertex
    edge = trichoplax.edge
    edge2 = trichoplax.edge2

    # cytoskeleton spring energy
    for j in 1:trichoplax.n_edges2[i]
        k = edge2[i,j]  # edge index of jth edge at ith vertex
        e = edge[k, :]     # vertex index of jth edge at ith vertex
        r = sqrt( ( v[e[1],1] - v[e[2],1] ) ^2 + ( v[e[1],2] - v[e[2],2] ) ^2 )
        Le = Le + r
        Es = Es + trichoplax.k2[]*(r - trichoplax.edgelength[k])^2
    end

    # cell turgor pressure
    Ep = 0.0
    tv = trichoplax.vertexcells
    tc = trichoplax.cell
    for j in 1:trichoplax.n_vertexcells[i]
        cellvertices = v[tc[tv[i,j],:],:]
        Ep = Ep + (cellvolume(cellvertices)-
                    trichoplax.volume[tv[i,j]])^2
    end


    # surface energy of external membranes
    # proportional to length of skin segments each side of skin
    Lx = 0.0
    if any(trichoplax.skin.==i)   # if i is a skin vertex
        j = findfirst(trichoplax.skin.==i)[]   # which skin vertex?
        snbr = trichoplax.skin_neighbour[j, :]
        Lx = sqrt( (v[i,1]-v[snbr[1],1])^2 + (v[i,2]-v[snbr[1],2])^2 +
                   (v[i,1]-v[snbr[2],1])^2 + (v[i,2]-v[snbr[2],2])^2)
    end

    # println(Es, ", ", trichoplax.σ[]*Lx)
     return (Es +  trichoplax.ρ[]*Ep + trichoplax.σ[]*Lx )
end


function localShapeEnergyGradient(dv::Float64, trichoplax::Trichoplax)

    v = trichoplax.vertex
    n = size(v)
    ∇E = fill(0.0, n )

    for i in 1:n[1]
        for j in 1:2
            v_save = v[i,j]
            v[i,j] = v[i,j] + dv  # + perturb (i,j)th coordinate
            ∂Eplus = localShapeEnergy(i, trichoplax)
            v[i,j] = v[i,j] - 2.0*dv  # - perturb (i,j)th coordinate
            ∂Eminus = localShapeEnergy(i, trichoplax)
            ∇E[i,j] = (∂Eplus-∂Eminus)/(2.0*dv) # ∂E/∂v_ij
            v[i,j] =v_save      # put (i,j)th coordinate back
        end
    end
    ∇E
end

function morph(trichoplax, rate::Float64 = .005, nsteps::Int64 = 25)

  @inbounds for t = 1:nsteps
      ∇ = localShapeEnergyGradient(1e-4*trichoplax.celldiameter, trichoplax)
      trichoplax.vertex[:,:] = trichoplax.vertex - rate*∇
  end

  #trichoplax = cellverticesfromskeleton(trichoplax)

  return trichoplax

end


#     ncells = size(cell,1)
#
#     #number of cell layers in body
#     nlayers = Int64((3+sqrt(9+12*(ncells-1)))/6)
#
#     #number of cells in map (outer celldepth layers)
#     nmap = ncells - (3*(nlayers-mapwidth)*(nlayers-mapwidth-1) + 1)
#
#     # duplicate body
#     mapcell = fill(0, nmap, 6)
#
#     nvertex = size(vertex,1)
#     @inbounds for i in 1:nmap  #size(cell,1)
#         @inbounds for j in 1:6
#             mapcell[i,j] = cell[ncells-i+1,j]
#         end
#     end
#
#     # minimum cell index in pCell
#     i0 = findmin(mapcell)[1]
#
#     # copy the required vertices, shift indices to coincide
#     #pVertex = copy(vertex[i0:end,:])
#     nmapvertex = nvertex-i0+1
#     mapvertex = fill(0.0, nmapvertex,2 )
#     mapcell = mapcell.-(i0-1)
#
#     # project vertices beyond skin
#     Rs = sqrt(vertex[end,1]^2 + vertex[end,2]^2) # skin radius
#     Rc = sqrt(vertex[2,1]^2 + vertex[2,2]^2)  # cell radius
#     @inbounds for i in 1:nmapvertex
#        R0 = sqrt(vertex[i0+i-1,1]^2 + vertex[i0+i-1,2]^2)
#        R1 = Rs + (Rs-R0)*(Rs/R0)^2
#        mapvertex[i,:] = vertex[i0+i-1,:].*R1/R0
#     end
#
#     (mapvertex, mapcell)
# end
#
# function makereceptivefields(vertex, cell, clayer, mapwidth)
#     # construct receptive field centers
#     # by reflecting centroids of mapping cells
#     # into annular region of width senseRange
#
#     ncells = size(cell,1)
#
#     #number of cell layers in body
#     nlayers = Int64((3+sqrt(9+12*(ncells-1)))/6)
#
#     # number of cells in map (outer celldepth layers)
#     nmap = ncells - (3*(nlayers-mapwidth)*(nlayers-mapwidth-1) + 1)
#
#     # index of cell (centroid) preceding first cell in map
#     i0 = 3*(nlayers-mapwidth)*(nlayers-mapwidth-1)+1
#
#     # array of rf centers
#     rfcenter = fill(0.0, nmap,2)
#
#     # project vertices beyond skin
#     Rs = sqrt(vertex[end,1]^2 + vertex[end,2]^2) # skin radius
#     Rc = sqrt(vertex[2,1]^2 + vertex[2,2]^2)  # cell radius
#     @inbounds for i in 1:nmap
#         x = mean(vertex[cell[i0+i,:],1])
#         y = mean(vertex[cell[i0+i,:],2])
#         R0 = sqrt(x^2 + y^2)
#         R1 = Rs + (Rs-R0)*(Rs/R0)^2
#        rfcenter[i,:] = [x*R1/R0 y*R1/R0]
#     end
#
#     return rfcenter
# end

# end  # module Placozoan
