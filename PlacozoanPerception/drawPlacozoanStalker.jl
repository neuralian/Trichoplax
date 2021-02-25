# draw placozoan
#

fig = Figure(resolution = (800,800))
        
ax = Axis(fig, backgroundcolor = "#918d8c")


matPlot = poly!(decompose(Point2f0, Circle(Point2f0(0,0), mat_radius)), 
   strokecolor = "#bac7b3", color = "#e1f2d8")

psPlot = poly!(decompose(Point2f0, Circle(Point2f0(0,0), prey.radius)), 
   strokecolor = "#bac7b3", color = "#ffedf3")
psGutPlot = poly!(decompose(Point2f0, Circle(Point2f0(0,0), prey.gutradius)), 
   strokecolor = "#c9b5c9", color = "#debfd4")

receptor_plt = scatter!(ax, prey.receptor.x, prey.receptor.y ,
   markersize=10.0,
   color=[(i==2) ? "#fccb28" : "#b3a7a4" for i in 1:prey.receptor.N],
   strokecolor=:black, strokewidth=0.25)

ppPlot = poly!(decompose(Point2f0, Circle(Point2f0(250.,250.), predator.radius)), 
   strokecolor = "#572e24", color = ("#855d52", .2))

Φ = 0.8175*π/2
mx = prey.mcell.d*cos(Φ)
my = prey.mcell.d*sin(Φ)
mc_coords = decompose(Point2f0, Circle(Point2f0(mx,my), prey.mcell.r))
mcellPlot = poly!(mc_coords, 
    strokecolor = "#a1320d", color = "#eb4710")


# project mc coords into the world
N = length(mc_coords)
rf_coords = fill(Point2f0(0,0), N)
for i in 1:N
    p = mc_coords[i]
    Φ = atan(p[2], p[1])
    r = sqrt(p[1]^2 + p[2]^2)
    r1 = (prey.radius - r)*(mat_radius -prey.radius)/prey.marginwidth + prey.radius
    rf_coords[i] = Point2f0(r1*cos(Φ), r1*sin(Φ))
end

rfPlot = poly!(rf_coords, 
    strokecolor = "#a1320d", color = "#eb4710")

fig[1,1] = ax
xlims!(ax, [0, 400])
ylims!(ax, [0, 400])
hidespines!(ax)
hidedecorations!(ax)

display(fig)