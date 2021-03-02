# draw placozoan
#

using CairoMakie, ColorSchemes, Colors
include("makieThemeBlue.jl") 

mat_radius = 400
r = 401.


figA = Figure(resolution = (801,801))

axA = Axis(figA, backgroundcolor = :white)

RF4 = copy(prey.receptor.pOpen[4])
for i in -prey.radius:prey.radius
   for j in -prey.radius:prey.radius
       r = sqrt(i^2+j^2)
       if (r<prey.radius) & (r > (prey.radius-prey.marginwidth))  # in marginal zone
         # project map location to real-world location
        # R = p.radius + (p.radius - r*(p.observer.maxRange - p.radius))/p.marginwidth 
 
         R = (prey.radius - r)*(prey.observer.maxRange-prey.radius)/ prey.marginwidth + prey.radius
                 
         #println(i, ", ", j, ", ", r, ", ", R, ", ", R/r)
         iproj = Int64(round(i*R/r))   
         jproj = Int64(round(j*R/r))
         # copy likelihood from world to map
         RF4[Int64(i),Int64(j)] = RF4[iproj,jproj]
       end
     end
   end



#heatmap!(OffsetArrays.no_offset_view(prey.receptor.pOpen[4]), colormap = :hot)
contour!(OffsetArrays.no_offset_view(RF4), color = :black, levels = [.105, .11, .125, .15, .2])

psPlot = poly!(decompose(Point2f0, Circle(Point2f0(r,r), prey.radius)), 
strokecolor = :red, color = RGBA(0,0,0,0), strokewidth = 2)

matPlot = poly!(decompose(Point2f0, Circle(Point2f0(r,r), mat_radius)), 
strokecolor = :green, color  = RGBA(0,0,0,0), strokewidth = 2)

# matPlot = poly!(decompose(Point2f0, Circle(Point2f0(r,r), mat_radius)), 
#    strokecolor = :white, color = RGBA(0,0,0,0))


# psGutPlot = poly!(decompose(Point2f0, Circle(Point2f0(r,r), prey.gutradius)), 
#    strokecolor = "#c19ba8", color = "#c19ba8")

receptor_plt = scatter!(axA, [r.+prey.receptor.x[4]], [r.+prey.receptor.y[4]] ,
   markersize=8.0, color = "#fccb28", 
  # color=[(i==4) ? "#fccb28" : "#b3a7a4" for i in 1:prey.receptor.N],
   strokecolor=:black, strokewidth=0.25)


figA[1,1] = axA
# xlims!(axA, [mat_radius 2*mat_radius])
# ylims!(axA,  [mat_radius 2*mat_radius])
hidespines!(axA)
hidedecorations!(axA)

display(figA)

save("LF_fig.png", figA, px_per_unit = 3 )

# figB = Figure(resolution = (801,801))

        
# axB = Axis(figB, backgroundcolor = "#918d8c")

# heatmap!(OffsetArrays.no_offset_view(1.0.-prey.receptor.pOpen[4]), colormap = :copper)


# matPlot = poly!(decompose(Point2f0, Circle(Point2f0(r,r), mat_radius)), 
#    strokecolor = :white, color = RGBA(0,0,0,0))

# psPlot = poly!(decompose(Point2f0, Circle(Point2f0(r,r), prey.radius)), 
#    strokecolor = "#e3bdb5", color = "#e3bdb5")
# psGutPlot = poly!(decompose(Point2f0, Circle(Point2f0(r,r), prey.gutradius)), 
#    strokecolor = "#c19ba8", color = "#c19ba8")

# receptor_plt = scatter!(axB, [r.+prey.receptor.x[4]], [r.+prey.receptor.y[4]] ,
#    markersize=10.0, color = "#1f59ba", 
#   # color=[(i==4) ? "#fccb28" : "#b3a7a4" for i in 1:prey.receptor.N],
#    strokecolor=:black, strokewidth=0.25)


# figA[1,1] = axB
# # xlims!(ax, [0, mat_radius])
# # ylims!(ax, [0, mat_radius])
# hidespines!(axB)
# hidedecorations!(axB)

# display(figB)