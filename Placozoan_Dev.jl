# Placozoan_Dev.jl

using Makie
using Colors
#using LinearAlgebra
using Statistics
using Placozoan


function draw(trichoplax::Trichoplax, color=:black, linewidth = 1.0)

    @inbounds for i in 1:size(trichoplax.cell,1)
        lines!(trichoplax.vertex[trichoplax.cell[i,[1:6; 1]],1],
               trichoplax.vertex[trichoplax.cell[i,[1:6; 1]],2],
                color = color, linewidth=linewidth)
    end

end

function drawmap(trichoplax, color = :blue, linewidth = 0.5)
  drawcells(trichoplax.mapvertex, trichoplax.mapcell, color, linewidth)

end

function drawskeleton(trichoplax::Trichoplax,
         color = RGB(.25,.65,.25), linewidth = 0.25)
  for link in 1:size(trichoplax.skeleton.edge,1)
      lines!(trichoplax.skeleton.vertex[trichoplax.skeleton.edge[link, :],1],
             trichoplax.skeleton.vertex[trichoplax.skeleton.edge[link, :],2],
             color=color, linewidth=linewidth)
  end
end


function paint(trichoplax, color = :red)

    nVertex = size(trichoplax.vertex, 1)
    nCell = size(trichoplax.cell, 1)
    vertex = fill(0.0, nVertex+nCell, 2) # will add centre vertex to each cell
    for i in 1:nVertex
        vertex[i,:] = trichoplax.vertex[i,:]
    end
    facet = fill(0,6*nCell,3) # 6 triangle facets per cell
    nFacet = 0
    for cell in 1:nCell
        vertex[nVertex+cell, :] = mean(vertex[trichoplax.cell[cell,:],:], dims=1)
        for i in 1:6
            nFacet = nFacet + 1
            i0 = trichoplax.cell[cell,i]
            i1 = trichoplax.cell[cell,i%6+1]
            facet[nFacet, :] = [nVertex+cell i0 i1]
        end
    end

    # render cells
    poly!(vertex, facet, color = rand(size(vertex,1)),
            colormap = :blues)


end


function skeletonEnergy(trichoplax::Trichoplax)
    #  potential energy in skeleton deformation + surface energy.
    #   skeleton edges are linear springs Es = (1/2)k(r-ro)^2
    #   surface energy σ per unit length of exposed membrane

# probably redundant??

    k2 = 1.0  # half of cytoskeleton spring constant (k/2)
    σ = 1.0  # surface energy density

    Es = 0.0
   # elastic energy in skeleton deformation
   for i in 1:size(trichoplax.skeleton.edge,1)
       r = sqrt(sum(
            ( trichoplax.skeleton.vertex[trichoplax.skeleton.edge[i,1],:] -
              trichoplax.skeleton.vertex[trichoplax.skeleton.edge[i,2],:] ).^2))
       Es = Es + k2*(r - trichoplax.skeleton.edgelength)^2
   end

   # surface energy of external membranes
   Lx = 0.0    # initialize external surface length to 0
   # x & y coords of external vertices
   x = mean(trichoplax.skeleton.vertex[trichoplax.skeleton.distalΔ[:,:], 1],
            dims = 2)
   y = mean(trichoplax.skeleton.vertex[trichoplax.skeleton.distalΔ[:,:], 2],
            dims = 2)
   # length of external surface, sum lengths of external edge segments
   n = length(x)
   for i in 2:n
       Lx = Lx + sqrt( (x[i]-x[i-1]).^2 + (y[i]-y[i-1]).^2)
   end
   # include edge from last to first vertex (close the loop)
   Lx = Lx + sqrt( (x[n]-x[1]).^2 + (y[n]-y[1]).^2)

   # total energy is elastic + surface energy
   println(Es, ", ", σ*Lx)
    return (Es +  σ*Lx)

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
    # x & y coords of external vertices
    x = mean(vertex[trichoplax.skeleton.distalΔ[:,:]], dims = 2)
    y = mean(vertex[trichoplax.skeleton.distalΔ[:,:].+nv], dims = 2)
    # length of external surface, sum lengths of external edge segments
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
        ∇[i] = skeletonEnergy(vertex, trichoplax) - E0
        vertex[i] = vertex[i] - dv
    end
    ∇
end


function morph(trichoplax)
  s = size(trichoplax.skeleton.vertex)
  for t = 1:5
      v = trichoplax.skeleton.vertex[:]
      ∇ = skeletonEnergyGradient(.1, v, trichoplax)
      trichoplax.skeleton.vertex[:,:] = reshape(v + ∇, s...)
  end

  trichoplax
end



# MAIN
bodylayers = 2 # number of body cell layers
mapdepth = 1     # map layers
celldiam = 10.0




@time trichoplax = Trichoplax(bodylayers, celldiam, mapdepth)
trichoplax.k2[] = 1.0   # cytoskeleton spring constant /2
trichoplax.σ[] = 1.0    # cell surface energy density

# Draw
s = Scene(resolution = (800,800), scale_plot = false)
draw(trichoplax, :red)
#drawmap(trichoplax, :orange, 1.0)
drawskeleton(trichoplax, RGB(.0, .65, .0), .5)

display(s)
