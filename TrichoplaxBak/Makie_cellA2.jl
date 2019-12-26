using Makie
using Colors
using Parameters
using LinearAlgebra
using BenchmarkTools

@with_kw mutable struct Cell

   # GEOMETRY
   # vertices, 3 hexagons (apical, medial, basal)
   # point in centre of apical and basal face
   v::Array{Float64,2} = vcat(
     [0 0 1],  # apex
      hcat([ cos.(2*π*i/6) for i in 0:5],
           [sin.(2*π*i/6) for i in 0:5], ones(6,1)), # apical face
      hcat([ cos.(2*π*i/6) for i in 0:5],
           [sin.(2*π*i/6) for i in 0:5], zeros(6,1)), # waist
      hcat([ cos.(2*π*i/6) for i in 0:5],
           [sin.(2*π*i/6) for i in 0:5], -ones(6,1)), # basal face
      [0 0 -1] # base
      )

   # face indices
   # each triplet defines triangular face, anticlockwise out.
   f::Array{UInt8,2} = UInt8[
       1 2 3 ; 1 3 4;  1 4 5; 1 5 6; 1 6 7; 1 7 2   # apical

       # lateral walls
       2 8 9 ; 2 9 3; 3 9 10; 3 10 4; 4 10 11; 4 11 5; 5 11 12;
       5 12 6; 6 12 13; 6 13 7; 7 13 8; 7 8 2;
       8 14 15; 8 15 9; 9 15 16; 9 16 10; 10 16 17; 10 17 11; 11 17 18;
       11 18 12; 12 18 19; 12 19 13; 13 19 14; 13 14 8;

       20 19 18; 20 18 17; 20 17 16; 20 16 15; 20 15 14; 20 14 19 # basal
       ]
  nfaces = size(f,1)



   # edges
   e = [    2 3; 3 4; 4 5; 5 6; 6 7; 7 2;
            2 8; 3 9; 4 10; 5 11; 6 12; 7 13;
            #8 9; 9 10; 10 11; 11 12; 12 13; 13 8;
            8 14; 9 15; 10 16; 11 17; 12 18; 13 19;
            14 15; 15 16; 16 17; 17 18; 18 19; 19 14
            ]

   edge = reshape(v[e',:], 2*size(e)[1], 3)
   facecolor = fill(RGBA(1.0, 0.0, 0.0, 0.5),nfaces)
   edgecolor = RGB(.4,0,0)
   edgewidth = 4;

   volume::Float64 = 0;

   # DYNAMICS
   surface_energy = 1e-6;  # approx surface energy of a cell, uJ.cm^-2

end


render(c::Cell) = begin
   # individual face coloring
   scene = mesh(c.v,c.f[1,:]', color = c.facecolor[1]);
   for iface = 2:size(c.f,1)
     mesh!(c.v, c.f[iface,:]', color = c.facecolor[iface]);
   end
   linesegments!(c.edge, color = c.edgecolor, linewidth = c.edgewidth);

  display(scene)

end

render(c::Cell, facecolor) = begin
   # all faces the same color, renders about 5x faster
   scene = mesh(c.v,c.f, color = facecolor);
   linesegments!(c.edge, color = c.edgecolor, linewidth = c.edgewidth);

  display(scene)

end

c= Cell()
@time render(c)
@time render(c, RGBA(1,0,0, 1))

cellvolume(c::Cell) = begin
   # volume of tetrahedron with 1 vertex at 0
   V = 0.0;
   for i in 1:size(c.f)[1]
      V += det(c.v[c.f[i,:],:])
   end
   return V/6.0
end


cellvolume(c)

#hello
