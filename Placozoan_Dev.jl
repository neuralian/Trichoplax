# Placozoan_Dev.jl

using Makie
using Colors
#using LinearAlgebra
using Statistics
using Placozoan


function draw(trichoplax::Trichoplax, color=:black, linewidth = 1.0)
  drawcells(trichoplax.vertex, trichoplax.cell, color, linewidth)
end

function drawmap(trichoplax, color = :blue, linewidth = 0.5)
  drawcells(trichoplax.mapvertex, trichoplax.mapcell, color, linewidth)

end

function drawskeleton(trichoplax::Trichoplax,
         color = RGB(.5,.5,.5), linewidth = 0.25)
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


function deformationEnergy(trichoplax)
    # mechanical energy stored in deformation of trichoplax skeleton
    # + surface energy.
    #   skeleton links are springs with rest length celldiam, Es = 1/2k(x-xo)^2
    #   surface energy of external cell membrane segment of length L, Ex = σL


    hk = 1.0  # half of cytoskeleton spring constant (k/2)
    σ = 100.0  # surface energy density

    Es = 0.0

    # aliases to simplify code
    vs = trichoplax.skeleton.vertex
    e = trichoplax.skeleton.edge
    r0 = trichoplax.celldiam

   # elastic energy in skeleton deformation
   for i in 1:size(e,1)
       r = sqrt(sum((vs[e[i,:],1] - vs[e[i,:],2]).^2))
       Es = Es + hk*(r - r0)^2
   end

   # surface energy of external membranes
   Ex = 0.0
   vt = trichoplax.vertex
   s = trichoplax.skinvertex
   nSkinVertices = length(s)
   j = sortperm(atan.(vt[s,1], vt[s,2])) # order of skin vertices anticlockwise
   for i in 2:nSkinVertices
       Ex = Ex + sqrt(sum((vt[s[j[i]],:] - vt[s[j[i-1]],:]).^2))
   end
   Ex = Ex + sqrt(sum((v[s[j[nSkinVertices]],:] - v[s[j[1]],:]).^2))
   Ex = σ*Ex

   E = Es + Ex

   return E

end

function deformationEnergy(trichoplax, vertex)
    # deformation energy as a function of skeleton and perimeter vertex coords



end


# MAIN
bodylayers = 3    # number of body cell layers
mapdepth = 1     # map layers
celldiam = 10.0


trichoplax = Trichoplax(bodylayers, celldiam, mapdepth)

# Draw
s = Scene(resolution = (800,800), scale_plot = false)
draw(trichoplax, :red)
#drawmap(trichoplax, :orange, 1.0)
drawskeleton(trichoplax, RGB(.5, .5, 1.0))

display(s)
