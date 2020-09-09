# PredPrey 03
# Predator and prey move; visualize in predator frame

using Makie

sceneWidth = 5000.0
sceneLimits = FRect(-1000.0, -1000.0, sceneWidth, sceneWidth)

predatorRadius = 500.   # body radius
preyRadius = 400.
c = 20.   # cell radius

preyLocation = (3000.0, 3000.0)

bacteriaDensity = 0.001  # bacteria /um^2

scene = Scene(resolution = (1000, 1000), limits = sceneLimits)

# time observable
t = Node(0.0)

nBacteria = sceneWidth*sceneWidth*bacteriaDensity

predator = poly!(scene,
       decompose(Point2f0, Circle(Point2f0(0.,0.), predatorRadius)),
       color = :gray, strokewidth = 1, strokecolor = :black)

prey = poly!(scene,
       lift(s->decompose(Point2f0, Circle(Point2f0(preyLocation),
       preyRadius)), t), color = :pink, strokewidth = 1, strokecolor = :black)

#bacteria = scatter!(scene, rand(nBacteria), rand(nBacteria))

N = 200
#for i in 1:N
record(scene, "test.mp4", 1:N) do i

    predatorStep = randn(2).*5.
    preyStep = randn(2).*10.
    global preyLocation = preyLocation .+ preyStep

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
