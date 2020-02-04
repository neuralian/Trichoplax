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
# using Colors

# export Trichoplax,
#                 discworld, neighbour, makebody, findperimetervertices,
#         smoothperimeter, makecellmap, makereceptivefields,
#         drawdelaunaydisc, drawcells, drawskeleton, findperimeteredges,
#         distalskeleton

struct Skeleton
    vertex::Array{Float64,2}          # ?x2 vertices
    #edge::Array{Int64,2}              # ?x2 edges forming Delaunay Δ-ation
    edgelength::Float64               # nominal edge length (= cell diameter)
    layercount::Array{Int64,1}        # vertices in ith layer
    neighbour::Array{Int64,2}         # 6 neighbours of each internal vertex
    # distalΔ::Array{Int64,2}         # ?x3 vertices of outer layer of triangles
end

struct Trichoplax
    nlayers::Int64
    skeleton::Skeleton
    vertex::Array{Float64}    # cell vertices
    cell::Array{Int64}        # index vertices for each cell
    triangle::Array{Int64,2}  # skeleton vertex triangles containing
                              #    cell vertices
    skinstart::Int64          # triangle[skinstart:end,:] gives skin vertices
    # edge::Array{Int64}                 # [i,:] index links between cells
    # skinvertex::Array{Int64}           #  skin vertex indices
    # mapdepth::Int64                    # number of cell layers in sensory map
    # mapvertex
    # mapcell
    # k2::Array{Float64,1}  # half of cytoskeleton spring constant (k/2)
    # σ::Array{Float64,1}   # surface energy density
    # nb fields are immutable
    #    Declaring k2 & σ as arrays allows access via e.g. trichoplax.k2[]
end

# normaldistribution = Normal()

numEdges(num_cells) = num_cells*6 + 9*num_cells*(num_cells-1)

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

function neighbours(v, layercount)
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

function plotneighbours(v, nbrs)
    for i in 1:size(nbrs,1)
        for j = 1:6
            plot!(vec([v[i,1] v[nbr[i,j],1]]), vec([v[i,2] v[nbr[i,j],2]]))
        end
    end

end

function makecells(skeleton)
    # cell vertices given skeleton vertices v & neighbours nbr
    # also returns index of vertices for each cell
    #  and index of skeleton triangles containing cell vertices
    #  so that cell vertices can be quickly updated when skeleton moves

    tol = .01*skeleton.edgelength

    nCells = size(skeleton.neighbour,1)
    nVertices = 6*(size(skeleton.layercount,1)-1)^2
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

    return (vertex, cell, triangle)

end

function updatecells(trichoplax)
    #    function updatecells(v,triangle)

    n = size(trichoplax.triangle, 1)
    vertex = fill(0.0, n, 2)
    for i in 1:n
        vertex[i,:] = sum(v[trichoplax.triangle[i,:],:], dims=1)/3.
    end
    return vertex
end


function skinvertex(trichoplax)
    # external vertices

    i0 = trichoplax.skinstart
    n = size(trichoplax.triangle, 1) - i0 + 1
    skinvertex = fill(0.0, n, 2)
    v = trichoplax.skeleton.vertex
    for i in 1:n
        skinvertex[i,:] = sum(v[trichoplax.triangle[i+i0-1,:],:], dims=1)/3.
    end
    return skinvertex
end



function draw(trichoplax::Trichoplax, color=:black, linewidth = .25)

    @inbounds for i in 1:size(trichoplax.cell,1)
        lines!(trichoplax.vertex[trichoplax.cell[i,[1:6; 1]],1],
               trichoplax.vertex[trichoplax.cell[i,[1:6; 1]],2],
                color = color, linewidth=linewidth)
    end

end


function skeletonEnergy(vertex::Vector{Float64}, trichoplax::Trichoplax)
#  potential energy in skeleton deformation + surface energy.
#   skeleton edges are linear springs Es = (1/2)k(r-ro)^2
#   surface energy σ per unit length of exposed membrane
#  skeleton energy with specified vertices
# (vertex = trichoplax.skeleton.vertex[:], i.e. nv x 2 flattened to 2*nv x 1
k2 = 1.0  # half of cytoskeleton spring constant (k/2)
σ = 1.0  # surface energy density
    Es = 0.0
    nv = size(trichoplax.skeleton.vertex,1)
    for i in 1:size(trichoplax.skeleton.edge,1)
        r = sqrt(  ( vertex[trichoplax.skeleton.edge[i,1]] -
                     vertex[trichoplax.skeleton.edge[i,2]]    ) ^2 +
                   ( vertex[trichoplax.skeleton.edge[i,1]+nv] -
                     vertex[trichoplax.skeleton.edge[i,2]+nv] ) ^2
                )
        Es = Es + k2*(r - trichoplax.skeleton.edgelength)^2

    end

    # surface energy of external membranes
    Lx = 0.0    # initialize external surface length to 0
    sv = skinvertex
    n = length(x)
    for i in 2:n
        Lx = Lx + sqrt( (x[i]-x[i-1]).^2 + (y[i]-y[i-1]).^2)
    end
    # include edge from last to first vertex (close the loop)
    Lx = Lx + sqrt( (x[n]-x[1]).^2 + (y[n]-y[1]).^2)

    #println(Es, ", ", σ*Lx)
     return (Es +  σ*Lx)


end

function skeletonEnergyGradient(dv, vertex::Vector{Float64}, trichoplax::Trichoplax)

    n = length(vertex)
    ∇ = fill(0.0, n )
    E0 = skeletonEnergy(vertex, trichoplax)
    for i in 1:n
        vertex[i] = vertex[i] + dv
        ∇[i] = (skeletonEnergy(vertex, trichoplax) - E0)/dv
        vertex[i] = vertex[i] - dv
    end
    ∇
end


function morph(trichoplax)
  s = size(trichoplax.skeleton.vertex)
  for t = 1:100
      v = trichoplax.skeleton.vertex[:]
      ∇ = skeletonEnergyGradient(.01, v, trichoplax)
      trichoplax.skeleton.vertex[:,:] = reshape(v - .005*∇, s...)
  end

  trichoplax
end


function Skeleton(n_cell_layers, cell_diameter)

    layercount = skeletonlayercount(n_cell_layers)
    vertex = skeletonvertices(layercount, cell_diameter)
    neighbour = neighbours(vertex, layercount)

    # distalΔ = distalskeleton(nbr, layer)

    return Skeleton(vertex, cell_diameter, layercount, neighbour)
end

# function findperimetervertices(cell, vlink, clayer)
#     # find vertices on edge of disc
#     n = length(clayer[end])
#     outerperimeter = fill(0, n+6) # linked only to perimeter vertices
#     cornerperimeter = fill(0, n)  # linked to 1 outer and 1 internal vertex
#     innerperimeter = fill(0, n)   # internal vertices linked to corner vertices
#
#
#     # outerperimeter (outer perimeter) vertices have no vlinks to internal vertices
#     no = 0
#     @inbounds for i in 1:n  # for each cell in outer layer
#         @inbounds for v in 1:6  # vertices with < 3  links ..
#             if sum(vlink[:]'.==cell[clayer[end][i],v])<3
#                 no = no + 1
#                 outerperimeter[no] = cell[clayer[end][i],v]
#             end
#         end
#     end
#
#     # cornerperimeter (inner perimeter) vertices
#     # are vlinked to outerperimeter vertices + 1 internal vertex
#     nc = 0
#     @inbounds for i in 1:length(outerperimeter)
#         @inbounds for j in 1:size(vlink,1) # for each vlink
#             if  (vlink[j,1]==outerperimeter[i])    &  # 1st of jth vlink links to ith vertex
#                  !any(outerperimeter.==vlink[j,2]) &   # and not to any other outerperimeter vertex
#                  !any(cornerperimeter.==vlink[j,2])     # and hasn't already been found
#                nc = nc+1
#                cornerperimeter[nc] = vlink[j,2]         # found a corner perimeter vertex
#            elseif (vlink[j,2]==outerperimeter[i]) &   # ditto for 2nd element in jth vlink
#                  !any(outerperimeter.==vlink[j,1]) &
#                  !any(cornerperimeter.==vlink[j,1])
#                nc = nc+1
#                cornerperimeter[nc] = vlink[j,1]
#             end
#         end
#     end
#
#     # internal vertices are vlinked to inner perimeter vertices
#     # (used to align lateral membranes of perimeter cells orthogonal to skin )
#     ni = 0
#     @inbounds for i in 1:length(cornerperimeter)
#         @inbounds for j in 1:size(vlink,1)
#             if (vlink[j,1]==cornerperimeter[i])    & # 1st in jth vlink is a corner
#                 !any(outerperimeter.==vlink[j,2]) & # AND 2nd is not a perimeter cell
#                 !any(innerperimeter.==vlink[j,2])   # AND not already found
#                ni = ni + 1
#                innerperimeter[ni] = vlink[j,2]    # found interior link to corner
#             elseif (vlink[j,2]==cornerperimeter[i])    & # 2nd in jth vlink is a corner
#                    !any(outerperimeter.==vlink[j,1]) & # AND 1st is not a perimeter cell
#                    !any(innerperimeter.==vlink[j,1])   # AND not already found
#                   ni = ni + 1
#                   innerperimeter[ni] = vlink[j,1]    # found interior link to corner
#             end
#         end
#     end
#     (outerperimeter, cornerperimeter, innerperimeter)
# end
#
# function smoothperimeter(vertex, op, cp, ip)
#     # move perimeter vertices onto a circle
#
#     No = length(op)
#     Nc = length(ip)
#     ro = fill(0.0, No)  # radii of outer vertices
#     ri = fill(0.0, Nc)  # radii of internal vertices
#
#     Ro = 0.0
#     @inbounds for i in 1:No
#         Ro  = Ro + sqrt(vertex[op[i],1]^2 + vertex[op[i],2]^2)
#         ro[i] = sqrt(vertex[op[i],1]^2 + vertex[op[i],2]^2)
#     end
#     Rc = 0.0
#     @inbounds for i in 1:Nc
#         Rc = Rc + sqrt(vertex[cp[i],1]^2 + vertex[cp[i],2]^2)
#         ri[i] = sqrt(vertex[ip[i],1]^2 + vertex[ip[i],2]^2)
#     end
#
#     # mean radius of perimeter vertices
#     R = (Ro/No + Rc/Nc)/2.0
#
#     # shift corner vertices to projection of inner vertices onto mean radius
#     @inbounds for i in 1:Nc
#         vertex[cp[i],1] = vertex[ip[i],1]*R/ri[i]
#         vertex[cp[i],2] = vertex[ip[i],2]*R/ri[i]
#     end
#
#     @inbounds for i in 1:No
#         vertex[op[i],1] = vertex[op[i],1]*R/ro[i]
#         vertex[op[i],2] = vertex[op[i],2]*R/ro[i]
#     end
#
#     vertex
# end

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
   (vertex, cell, triangle) = makecells(skeleton)

  # # perimeter vertices in 3 classes
  # perimetervertex = []#findperimetervertices(cell, link, layer)
  # # ordered vertices on perimeter
  # skinvertex = []#sort(vcat(perimetervertex[1],perimetervertex[2] ))
  #
  # (mapvertex,mapcell) = ([], []) #makecellmap(vertex, cell, layer, mapdepth);

  skinstart = 6*(n_cell_layers-1)^2+1 # index to first skin triangle

  return Trichoplax(n_cell_layers, skeleton, vertex, cell, triangle, skinstart)
end

# function makecellmap(vertex, cell, clayer, mapwidth)
#     # reflect outer layers of cells into environment
#     # by reflecting mapWidth layers of cells around body perimeter
#     # into annular region of width senseRange
#
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
