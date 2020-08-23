# PredPrey 03
# Predator and prey move; visualize in predator frame

using Makie
using Colors

sceneWidth  = 5000.0
sceneOffset = sceneWidth / 2.0
matRadius2  = (sceneWidth / 2.0)^2
sceneLimits = FRect(-sceneOffset, -sceneOffset, sceneWidth, sceneWidth)
x = LinRange(-sceneOffset, sceneWidth-sceneOffset, 1000)
y = LinRange(-sceneOffset, sceneWidth-sceneOffset, 1000)
# z = exp(-(ones(1000,1)*((x .- 2000.).^2)'  +  ((y.-2000.).^2)*ones(1,1000)) ./ 1.0e8)


predatorRadius = 500.   # body radius
predatorFieldRadius = 1000.
preyRadius = 400.
preyMargin = 100.
nReceptor = 24
c = 20.   # cell radius

preyLocation = (0.0, 0.0)
predatorLocation = (2000.0, 2000.0)



scene = Scene(resolution = (1000, 1000), limits = sceneLimits, show_axis=false)
# image!(x,y,z, colormap = :reds, alpha = 0.5)
# time observable
t = Node(0.0)

# scatter bacteria over a disc centred on the prey
# disc extends off initial scene so that more bacteria come into scene
# as the prey moves (scene moves in prey frame)
bacteriaRadius = sceneWidth*3.
bacteriaDensity = 5e-5 # bacteria /um^2
nBacteria = Int(round(bacteriaDensity*π*bacteriaRadius^2))
bacteriaPos = bacteriaRadius*(rand(nBacteria,2) .- 0.5)
bacteriaSize = 5.0*rand(nBacteria)

mat = poly!(scene,
       decompose(Point2f0, Circle(Point2f0(0.,0.), sqrt(matRadius2))),
       color = RGBA(.42, .66, .14, 1.0), strokewidth = 0, strokecolor = :black)

# plot bacteria
bacteria = scatter!(scene,
             lift(s->(bacteriaPos[:,1], bacteriaPos[:,2]), t),
             color = :teal, markersize = bacteriaSize,
             strokewidth = 0, strokecolor = :green)[end]

# Prey
prey = poly!(scene,
       decompose(Point2f0, Circle(Point2f0(0.,0.), preyRadius)),
       color = RGBA(0.9, 0.75, 0.65, 1.0), facealpha = 0.1,
       strokewidth = .5, strokecolor = RGBA(0.8, 0.65, 0.55, 1.0))
preyGut = poly!(scene,
      decompose(Point2f0, Circle(Point2f0(0.,0.), preyRadius-preyMargin)),
      color = RGBA(1., 0.75, 0.75, 1.0), facealpha = 0.1,
      strokewidth = 0.0, strokecolor = RGBA(0.9, 0.75, 0.65, 1.0))
receptorOffColor = [RGB(0.85, 0.75, 0.45) for i in 1:nReceptor]
receptorColor = [RGB(0.85, 0.75, 0.45) for i in 1:nReceptor]
receptor = scatter!(scene,
           [preyRadius.*(cos(2π*i/nReceptor), sin(2π*i/nReceptor))
               for i in 1:nReceptor],
            markersize = 8, color = receptorColor, strokewidth =0)[end]

preyStep = [0.0, 0.0]

predatorField = poly!(scene,
  lift(s->decompose(Point2f0, Circle(Point2f0(predatorLocation),
  predatorFieldRadius)), t),
  color = RGBA(.3, 1.0, .5, .5), strokewidth = 0, strokecolor = :black)
predator = poly!(scene,
       lift(s->decompose(Point2f0, Circle(Point2f0(predatorLocation),
       predatorRadius)), t),
       color = RGB(.6, 0.6, .5), strokewidth = .5, strokecolor = :black)

predatorStep = [0.0, 0.0]



# function moveBacteria(bacteria, dx, dy)
#
#     bacteria

N = 500
#for i in 1:N
record(scene, "test.mp4", 1:N) do i

    global predatorStep = 0.95*predatorStep +
                          0.05*randn(2).*50.0  .- predatorLocation ./ 10000.
    global preyStep = 0.95*preyStep + 0.05*randn(2).*50.
    global preyLocation = preyLocation .+ preyStep
    global predatorLocation = predatorLocation .+ predatorStep .+ preyStep
    global bacteriaPos += ones(nBacteria,1)*preyStep'
    bacteria.markersize =
         bacteriaSize.*(sum(bacteriaPos.^2,dims=2) .< matRadius2)[:]
    receptorColor = receptorOffColor[:]
    b = [Int(rand()[] > 0.5) for i in 1:nReceptor]
    receptorColor[findall(x->x==1, b)] .= RGB(1.0, 1.0, 0.0)
    receptor.color[] = receptorColor
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
