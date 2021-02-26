# draw placozoan
#

using CairoMakie, ColorSchemes, Colors
include("makieThemeBlue.jl") 
mat_radius = 300

fig = Figure(resolution = (800,800))
        
ax = Axis(fig, backgroundcolor = "#918d8c")


matPlot = poly!(decompose(Point2f0, Circle(Point2f0(0,0), mat_radius)), 
   strokecolor = "#808080", color = " #7c8557")

psPlot = poly!(decompose(Point2f0, Circle(Point2f0(0,0), prey.radius)), 
   strokecolor = "#bac7b3", color = "#f8ddb8")
psGutPlot = poly!(decompose(Point2f0, Circle(Point2f0(0,0), prey.gutradius)), 
   strokecolor = "#c9b5c9", color = "#f6cebf")



d_pred = 250.
ppPlot = poly!(decompose(Point2f0, Circle(Point2f0(d_pred,d_pred), predator.radius)), 
   strokecolor = "#572e24", color = ("#855d52", .2))
# separation line


receptor_plt = scatter!(ax, prey.receptor.x, prey.receptor.y ,
   markersize=10.0, color = "#b3a7a4",
  # color=[(i==2) ? "#fccb28" : "#b3a7a4" for i in 1:prey.receptor.N],
   strokecolor=:black, strokewidth=0.25)

   # list Φ and d for each m-cell to draw
d1 = prey.radius
d0 = prey.gutradius
Δ = d1-d0
M = 6
Φ = [π/32 π/10 π/6 π/4 0.34*π .44*π]
d = [d1-prey.mcell.r prey.mcell.d d0+3*Δ/4 d0+Δ/2+1.25 d0+Δ/4 d0+prey.mcell.r]
N = 64 # default number of points in circle decomposition
mcell_pts = [[ Point2f0(0,0) for _ in 1:N] for _ in 1:M]
rf_pts    = [[ Point2f0(0,0) for _ in 1:N] for _ in 1:M]
for i in 1:M

   mx = d[i]*cos(Φ[i])
   my = d[i]*sin(Φ[i])

   mcell_pts[i] = decompose(Point2f0, Circle( Point2f0(mx,my), prey.mcell.r) )
   poly!(mcell_pts[i],  strokecolor = "#a1320d", color = "#eb4710")

   # project mc coords into the world
   for j in 1:N
      p = mcell_pts[i][j]
      Ω = atan(p[2], p[1])
      r = sqrt(p[1]^2 + p[2]^2)
      r1 = (prey.radius - r)*(mat_radius -prey.radius)/prey.marginwidth + prey.radius
      rf_pts[i][j] = Point2f0(r1*cos(Ω), r1*sin(Ω))
   end
   poly!(rf_pts[i], strokecolor = "#a1320d", color = "#c4577c")
end

d3 = 205.
lines!([d1*cos(π/4), d3*cos(π/4)], [ d1*sin(π/4), d3*sin(π/4)], color =:red)
scatter!([ d3*cos(π/4)], [ d3*sin(π/4)], color = :red)

# scalebar
lines!([200; 250], [25; 25], color = :black, linewidth = 10)

fig[1,1] = ax
xlims!(ax, [0, mat_radius])
ylims!(ax, [0, mat_radius])
hidespines!(ax)
hidedecorations!(ax)

display(fig)

save("RF_fig.png", fig, px_per_unit = 3 )