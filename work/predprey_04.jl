# PredPrey 03
# Predator and prey move; visualize in predator frame

using Makie

sceneWidth  = 5000.0
sceneOffset = 1000.0
sceneLimits = FRect(-sceneOffset, -sceneOffset, sceneWidth, sceneWidth)

predatorRadius = 500.   # body radius
preyRadius = 400.
c = 20.   # cell radius

preyLocation = (0.0, 0.0)
predatorLocation = (3000.0, 3000.0)



scene = Scene(resolution = (1000, 1000), limits = sceneLimits)

# time observable
t = Node(0.0)

# scatter bacteria over a disc centred on the prey
# disc extends off initial scene so that more bacteria come into scene
# as the prey moves (scene moves in prey frame)
bacteriaRadius = sceneWidth*3.
bacteriaDensity = 1e-5 # bacteria /um^2
nBacteria = Int(round(bacteriaDensity*π*bacteriaRadius^2))
bacteriaPos = bacteriaRadius*(rand(nBacteria,2) .- 0.5)

prey = poly!(scene,
       decompose(Point2f0, Circle(Point2f0(0.,0.), preyRadius)),
       color = :pink, strokewidth = 1, strokecolor = :black)
preyStep = [0.0, 0.0]

predator = poly!(scene,
       lift(s->decompose(Point2f0, Circle(Point2f0(predatorLocation),
       predatorRadius)), t),
       color = :gray, strokewidth = 1, strokecolor = :black)

# plot bacteria
bacteria = scatter!(scene,
              lift(s->(bacteriaPos[:,1], bacteriaPos[:,2]), t),
              color = :green, markersize = 5)[end]

# function moveBacteria(bacteria, dx, dy)
#
#     bacteria

N = 100
#for i in 1:N
record(scene, "test.mp4", 1:N) do i

    predatorStep = randn(2).*5.
    global preyStep = 0.95*preyStep + 0.05*randn(2).*100.
    global preyLocation = preyLocation .+ preyStep
    global predatorLocation = predatorLocation .+ predatorStep
    global bacteriaPos += ones(nBacteria,1)*preyStep'

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
