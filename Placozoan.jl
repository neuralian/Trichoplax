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

export delaunaydisc, neighbours, makebody, findperimeter,
        smoothperimeter, makecellmap, makereceptivefields,
        drawdelaunaydisc, drawcells, drawskeleton


normaldistribution = Normal()

function delaunaydisc(nlayers, layerwidth)
    # returns (nucleus, dlink, layer)
    # nucleus = nCells x 2 coordinates of ith cell nucleus
    # dlink = ndlink x indices of cell nuclei defining the dlink
    # layer = ith row indexes cell nuclei in each layer
    # MGP Dec 2019

poswobble = .0025   # noise on cell location (re. layerwidth)
ncells = 3*nlayers*(nlayers-1)+1
nucleus = fill(0.0,ncells , 2)
dlink = fill(0, 6*ncells, 2)
layer = fill(Int64[],nlayers)
nlayer = vcat([1], [6*i for i in 1:nlayers-1]) # n cells in layer
nlinkpercell = fill(0, ncells, 1)  # number of Delaunay links per cell
icell = 1
nucleus[icell,:] = [0.0 0.0]  # cell at origin
layer[1,1] = [1]
idlink =0
w = 0
@inbounds for ilayer in 2:nlayers
    layer[ilayer,:] = [fill(0,nlayer[ilayer])]
    w = w + 2.0*(rand(1)[1] - 0.5)/nlayer[ilayer]
    @inbounds for j in 1:nlayer[ilayer]
        dw = rand(normaldistribution)
        while abs(dw)>3.0
            dw = rand(normaldistribution)
        end
        dw= .1*dw/(ilayer*layerwidth)
        dr = rand(normaldistribution)
        while abs(dr)>3.0
            dr = rand(normaldistribution)
        end
        dr = poswobble*dr
        icell = icell + 1
        layer[ilayer][j] = icell
        nucleus[icell,:] =
        (ilayer-1)*(1.0+dr)*layerwidth.*
        [ cos(2π*((j-1)/nlayer[ilayer])+w+dw)
          sin(2π*((j-1)/nlayer[ilayer])+w+dw) ]

        # Delaunay triangulation
        if j>1
            idlink = idlink + 1
            dlink[idlink, :] = [icell-1 icell]
            nlinkpercell[icell-1] = nlinkpercell[icell-1]+1
            nlinkpercell[icell] = nlinkpercell[icell]+1
            #lines!(nucleus[dlink[idlink,:],1], nucleus[dlink[idlink,:],2])
            #display(s)
        end
    end
    idlink = idlink + 1
    dlink[idlink, :] = [icell layer[ilayer][1]]
    nlinkpercell[icell] = nlinkpercell[icell]+1
    nlinkpercell[layer[ilayer][1]] = nlinkpercell[layer[ilayer][1]]+1
    #lines!(nucleus[dlink[idlink,:],1], nucleus[dlink[idlink,:],2])
    #display(s)
end
# print((ncells, icell))
# ncells = icell
# nucleus = nucleus[1:ncells,:]

# scatter!(nucleus[:,1], nucleus[:,2], markersize = layerwidth/4.)
# @inbounds for i in 1:nlayers
#      scatter!([nucleus[layer[i][1], 1]], [nucleus[layer[i][1],2]],
#      markersize = .5, color = :red);
#  end
# display(s)

# complete Delaunay triangles
@inbounds for j in 1:6      # center cell special case, 6 links
    idlink = idlink + 1
    dlink[idlink,:] = [1 j+1]
    nlinkpercell[1] = nlinkpercell[1]+1
    nlinkpercell[j+1] = nlinkpercell[j+1]+1
    #lines!(nucleus[dlink[idlink,:],1], nucleus[dlink[idlink,:],2])
    #display(s)
end

@inbounds for ilayer in 2:nlayers-1
    ilink = [ nlayer[ilayer+1] 1 2]
    link = layer[ilayer+1][ilink[1:3]]
    @inbounds for i in 1:length(layer[ilayer])
        #println((ilayer, link))
        @inbounds for j in 1:length(link)            #global idlink = idlink + 1
            idlink = idlink + 1
            dlink[idlink, :] = [layer[ilayer][i] link[j]]
            nlinkpercell[layer[ilayer][i]] = nlinkpercell[layer[ilayer][i]]+1
            nlinkpercell[link[j]] = nlinkpercell[link[j]]+1
            # lines!(nucleus[dlink[idlink,:],1], nucleus[dlink[idlink,:],2])
            # display(s)
            # sleep(.05)
        end
        #println((link[end],nlinkpercell[layer[ilayer][i+1]] ))
        if i<length(layer[ilayer])
            link = link[end] .+ collect(0:(5-nlinkpercell[layer[ilayer][i+1]]))
        end
        #println(Link)
    end
end
dlink = dlink[1:idlink,:]

return (nucleus, dlink, layer)

end


function neighbours(nucleus, dlink, layer)
  # index the six neighbour nuclei of each nucleus (excepting outer layer)
  # given output from delaunaydisc()...   (nb splatted argument)
  # outer nucleus layer
  # MGP Dec 2019

  ncells = size(nucleus, 1)
  neighbour = fill(0, ncells, 6)   # neighbour indices for each cell
  n_neighbour = fill(0, ncells)       # number of neighbours found

  # find neighbours of each cell
  @inbounds for i in 1:size(dlink,1)
    a = dlink[i,1]
    b = dlink[i,2]
    n_neighbour[a] = n_neighbour[a]+1
    n_neighbour[b] = n_neighbour[b]+1
    neighbour[a, n_neighbour[a]] = b
    neighbour[b, n_neighbour[b]] = a
  end
  # put neighbours in counterclockwise order
  # nb row 1 is already ccw, skip external nuclei (less than 6 links)
  @inbounds for i in 2:(ncells-length(layer[end]))
    dx = nucleus[neighbour[i,:],1].-nucleus[i,1]
    dy = nucleus[neighbour[i,:],2].-nucleus[i,2]
    order = sortperm(atan.(dy,dx))
    neighbour[i,:] = neighbour[i, order]
  end
  neighbour
end

function makebody(nucleus, neighbour, dlayer)
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


  ncells = dlayer[end-1][end]       #
  vertex = fill(0.0, 6*ncells, 2 )  # vertex coordinates
  cell = fill(0, ncells, 6)         # vertex indices for each cell
  vlink = fill(0, 6*ncells, 2)       # vertex indices for each vlink

  nvertex = 0
  nvlink = 0
  @inbounds for i in 1:ncells
      for j in 1:6
        k = neighbour[i,j]
        m = neighbour[i, j%6+1]
        x = (nucleus[i,1] + nucleus[k,1] + nucleus[m,1])/3.
        y = (nucleus[i,2] + nucleus[k,2] + nucleus[m,2])/3.
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

        # build list of vlinks (cell membrane segments)
        vlink_exists = false
        if j>1
            e = sort([cell[i,j-1] cell[i,j]], dims=2)  # vlink
            @inbounds for i1 in 1:nvlink
                if vec(e)==vlink[i1,:]
                    vlink_exists = true
                    break
                end
            end # for i1
            if !vlink_exists
            nvlink = nvlink + 1
            vlink[nvlink,:] = e
            end # !vlink_exists
        end # if j>1
        # 6th vlink links last and first vertex
        vlink_exists = false
        if j==6
            e = sort([cell[i,j] cell[i,1]], dims=2)
            for i1 in 1:nvlink
                if vec(e)==vlink[i1,:]
                    vlink_exists = true
                    break
                end
            end # for i1
            if !vlink_exists
                nvlink = nvlink + 1
                vlink[nvlink,:] = e
            end # !vlink_exists
        end # j==6
    end # for j
    end # for i
  (vertex[1:nvertex,:], cell, vlink[1:nvlink, :], dlayer[1:(end-1)])
end #cells()


function findperimeter(cell, vlink, clayer)
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


function drawdelaunaydisc(dlink,color = RGB(0.25, 0.25, 0.25), linewidth = 0.5)
    @inbounds for i = 1:size(dlink, 1)
      lines!(
        cellcentroid[dlink[i, :], 1],
        cellcentroid[dlink[i, :], 2],
        color = color, linewidth = linewidth
        )
    end
end

function drawcells(vertex, cell, color = RGB(.45, .25, 0.0), linewidth=1.0)
    @inbounds for i in 1:size(cell,1)
        lines!(vertex[cell[i,[1:6; 1]],1], vertex[cell[i,[1:6; 1]],2],
                color = color, linewidth=linewidth)
    end
end



function drawskeleton(vertex, vlink)
    # draw links between cell vertices (cytoskeleton)
    for i in 1:size(vlink, 1)
        lines!(vertex[vlink[i,:],1], vertex[vlink[i,:],2])
    end
end



end  # module Placozoan
