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
#z = exp(-(ones(1000,1)*((x .- 2000.).^2)'  +  ((y.-2000.).^2)*ones(1,1000)) ./ 1.0e9)


predatorRadius = 500.   # body radius
predatorFieldRadius = 1000.
preyRadius = 400.
preyMargin = 100.
nReceptor = 24
c = 20.   # cell radius

preyLocation = (0.0, 0.0)
predatorLocation = (2000.0, 2000.0)

pOpenState = fill(0.0, nReceptor, length(x), length(y))

# Johnson-Nyquist noise
kB = 1.38e-23
T = 300.
Ω = 20.e6   # 20MΩ
Δf = 1.0e3   # 1KHz
σ = sqrt(4.0*kB*T*Δf)

# channel open probability
kT = kB*T
v0 = -σ*log(0.1/(1.0-0.1))
p(e) =  1.0./(1 .+ exp.(-(e.-v0)/σ))

# function computes field strength at radius r from dipole source
# μV/cm
E_(r) = 1.0e6*2π*ρ*I*δ./r.^3
# array holds computed field strength from edge of Trichoplax
E = fill(0.0, nD)

# compute field strength at distance d from edge of body
function fieldStrength(sourceRadius, sourceMargin, nD, E)
  for a in cellDiam:cellDiam:(sourceRadius-cellDiam)
    n = round(2π*a/cellDiam)
    x = [ a*cos(2π*i/n) for i in 1:n]
    y = [ a*sin(2π*i/n) for i in 1:n]
    for d in 1:nD
      r = sqrt.(((d.+sourceRadius.-x).^2 + y.^2)).*1.0e-4
      E[d] = E[d] + sum(E_(r))
    end
    #scatter!(x, y, markersize = cellDiam/16)
  end
  return E
end

E = fieldStrength(preyRadius, preyMargin, nD, E)

function voltageMicrovolts(d, E)
  V = cumsum(E)*1.0e-4
  return V[end].-V[Int.(round.(d))]
end

V = voltageMicrovolts(1:nD, E)

function pOpen(d, V)
   i = Int(round(d)) + 1
   if i > length(V)
     i = length(V)
   end
   return p(V[i]*1.0e-6)
 end

# open state probability for source edge at (x,y) receptor at (x0,y0)
function pOpen(iReceptor, x0, y0)
  for i in 1:length(x)
    for j in 1:length(y)
       pOpenState[iReceptor, i,j] = pOpen(sqrt((x[i]-x0)^2 + (y[j]-y0)^2), V)
     end
   end
 end

scene = Scene(resolution = (1000, 1000), limits = sceneLimits, show_axis=false)
# surface!(scene,x,y,log.(abs.(z)),
#    color = r[Int.(round.(log.(abs.(z))))], overdraw = true, transparency=true)


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
receptorLocation = [preyRadius.*(cos(2π*i/nReceptor), sin(2π*i/nReceptor))
    for i in 1:nReceptor]
receptor = scatter!(scene, receptorLocation ,
            markersize = 8, color = receptorColor, strokewidth =0)[end]

preyStep = [0.0, 0.0]

# predatorField = poly!(scene,
#   lift(s->decompose(Point2f0, Circle(Point2f0(predatorLocation),
#   predatorFieldRadius)), t),
#   color = RGBA(.3, 1.0, .5, .5), strokewidth = 0, strokecolor = :black)
predator = poly!(scene,
       lift(s->decompose(Point2f0, Circle(Point2f0(predatorLocation),
       predatorRadius)), t),
       color = RGB(.6, 0.6, .5), strokewidth = .5, strokecolor = :black)

predatorStep = [0.0, 0.0]

for iReceptor in [1, 2, 24 ]
  pOpen(iReceptor, receptor[1][][iReceptor][1], receptor[1][][iReceptor][2])
end

# function moveBacteria(bacteria, dx, dy)
#
#     bacteria

N = 50
#for i in 1:N
b = fill(0.0, nReceptor)
record(scene, "test.mp4", 1:N) do i

    global predatorStep = 0.95*predatorStep +
                          0.05*randn(2).*50.0  .- predatorLocation ./ 5000.
    global preyStep = 0.95*preyStep + 0.05*randn(2).*50.
    global preyLocation = preyLocation .+ preyStep
    global predatorLocation = predatorLocation .+ predatorStep .+ preyStep
    global bacteriaPos += ones(nBacteria,1)*preyStep'
    bacteria.markersize =
         bacteriaSize.*(sum(bacteriaPos.^2,dims=2) .< matRadius2)[:]
    receptorColor = receptorOffColor[:]

    for j = 1:nReceptor
       range = sqrt( (predatorLocation[1] -receptorLocation[j][1])^2  +
                    (predatorLocation[2] -receptorLocation[j][2])^2 ) -
                     predatorRadius
       b[j] = Int(rand()[] < pOpen(range, V))
    end

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
q = pOpenState[1, :,:].*(1.0 .-pOpenState[2, :,:]).*(1.0 .-pOpenState[24, :,:])
contour!(x,y,q, colormap = :reds, levels = maximum(q)*[.5 , .9])
