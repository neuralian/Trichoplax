using Makie


R = 500.   # body radius
c = 20.   # cell radius

scene = Scene(resolution = (1000, 1000), limits = FRect(-5*R, -5*R, 10*R, 10*R))

# time observable
t = Node(0.0)

predatorLocation = (1000.0, 1000.0)

prey = poly!(scene,
       decompose(Point2f0, Circle(Point2f0(-1000.,-1000.), R)),
       color = :pink, strokewidth = 1, strokecolor = :black)

predator = poly!(scene,
       lift(s->decompose(Point2f0, Circle(Point2f0(predatorLocation),
       1.25*R)), t),
       color = :gray, strokewidth = 1, strokecolor = :black)

N = 100
#for i in 1:N
record(scene, "test.mp4", 1:N) do i

    predatorStep = randn(2).*5.
    global predatorLocation = predatorLocation .+ predatorStep

    t[] = 2*π*i/N
    #sleep(0.001)
    #display(scene)
end


# for r in c:c:(R-c)
#
#   n = round(2π*r/c)
#   ϕ = rand()
#   θ = 2π*((1:n)/n .+ ϕ)
#   scatter!((r)*sin.(θ), (r)*cos.(θ), markersize = c/8)
#
# end

display(scene)
