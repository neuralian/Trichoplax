using Makie


R = 500.   # body radius
c = 20.   # cell radius

scene = Scene(resolution = (1000, 1000), limits = FRect(-5*R, -5*R, 10*R, 10*R))

body = poly!(scene, decompose(Point2f0, Circle(Point2f0(0,0), R)),
            color = :lightgray, strokewidth = 1, strokecolor = :black)

# for r in c:c:(R-c)
#
#   n = round(2π*r/c)
#   ϕ = rand()
#   θ = 2π*((1:n)/n .+ ϕ)
#   scatter!((r)*sin.(θ), (r)*cos.(θ), markersize = c/8)
#
# end

display(scene)
