# PredPrey 03
# Predator and prey move; visualize in predator frame

using Makie
using Colors
using PerceptualColourMaps

sceneWidth  = 5000.0
sceneOffset = sceneWidth / 2.0
matRadius2  = (sceneWidth / 2.0)^2
sceneLimits = FRect(-sceneOffset, -sceneOffset, sceneWidth, sceneWidth)
sceneResolution = 1001  # nb must be odd, to ensure grid point @ centre
x = LinRange(-sceneOffset, sceneWidth-sceneOffset, 1001)
y = LinRange(-sceneOffset, sceneWidth-sceneOffset, 1001)
Nx = length(x)
Ny = length(y)
#z = exp(-(ones(1000,1)*((x .- 2000.).^2)'  +  ((y.-2000.).^2)*ones(1,1000)) ./ 1.0e9)


predatorRadius = 500.   # body radius
predatorMargin = 100.
preyRadius = 400.
preyMargin = 100.
nReceptor = 16  # multiple of 4
maxRange = 5.0e3 # max range (from centre) in μm
cellDiam = 20.   # cell radius

nD = Int(maxRange-preyRadius)

#preyLocation = (0.0, 0.0)
predatorLocation = (2000.0, 2000.0)

LikelihoodLookup = fill(1.0e-10, nReceptor, Nx, Ny)
LikelihoodArray = fill(0.0, Nx, Ny)

# indices of origin (0,0) point in LikelihoodArray
i0 = Int(ceil(length(x)/2.0))
j0 = Int(ceil(length(y)/2.0))

# source
ρ = 25.0   # Resisitivity of seawater 25Ω.cm
δ = 20.e-6*100.  # dipole separation 10μm in cm
I = 0.1e-12     # dipole current 0.1pA

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

# receptorState[i] == 1 if receptor i is active
receptorState = fill(0, nReceptor)

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

E = fieldStrength(predatorRadius, predatorMargin, nD, E)

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
       LikelihoodLookup[iReceptor, i,j] = pOpen(sqrt((x[i]-x0)^2 + (y[j]-y0)^2), V)
     end
   end
 end

 # function pOpen(iReceptor, x0, y0)
 #   for i in 1:length(x)
 #     for j in 1:length(y)
 #        LikelihoodLookup[iReceptor, i,j] = pOpen(sqrt((x[i]-x0)^2 + (y[j]-y0)^2), V)
 #      end
 #    end
 #  end

 function likelihood()

    LikelihoodArray .= 1.0
    for i = 1:nReceptor
      if receptorState[i]==1
        LikelihoodArray .*= LikelihoodLookup[i,:,:]
      else
        LikelihoodArray .*= (1.0 .- LikelihoodLookup[i,:,:])
      end
    end

   xx = findall(abs.(x).<=preyRadius)
   yy = findall(abs.(y).<=preyRadius)


    for i in xx
      for j in yy
        if x[i]^2 + y[j]^2 < preyRadius^2
            LikelihoodArray[i, j] = 0.0
          end
      end
    end

    LikelihoodArray ./= maximum(LikelihoodArray)

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

predator = poly!(scene,
      lift(s->decompose(Point2f0, Circle(Point2f0(predatorLocation),
      predatorRadius)), t),
      color = RGBA(.6, 0.6, .5, .75), strokewidth = .25, strokecolor = :black)

predatorStep = [0.0, 0.0]

# likelihood contours
receptorLocation = [preyRadius.*(cos(2π*i/nReceptor), sin(2π*i/nReceptor))
    for i in 1:nReceptor]
for j in 1:nReceptor
  receptorState[j] = Int(rand()[] < 0.5 )
end


# for each receptor construct lookup table
# for likelihood of predator at (x,y) given receptor state.
# Computed for receptors in first quadrant only,
# likelihoods for other quadrants are obtained by rotating this 90deg x 3
Nq = nReceptor ÷ 4  # number of receptors per quadrant
for iReceptor in 1:Nq
  pOpen(iReceptor, receptorLocation[iReceptor][1], receptorLocation[iReceptor][2])
end
# 2nd-4th quadrants
for iReceptor in 1:Nq
  # 2nd
  for i in 1:Nx
    for j in 1:Ny
      LikelihoodLookup[Nq+iReceptor, i,j] =
                            LikelihoodLookup[iReceptor, j,Nx+1-i]

    end
  end
  # 3rd
  for i in 1:Nx
    for j in 1:Ny
      LikelihoodLookup[2*Nq+iReceptor, i,j] =
                           LikelihoodLookup[iReceptor, Nx+1-i,Ny+1-j]
    end
  end
  # 4th
  for i in 1:Nx
    for j in 1:Ny
      LikelihoodLookup[3*Nq+iReceptor, i,j] =
                            LikelihoodLookup[iReceptor, Ny + 1 - j,i]
    end
  end
end


likelihood()
LhdPlot = contour!(scene, x,y,lift(s->LikelihoodArray, t),
         colormap = :Reds,
         levels = [0.1, .5 , .9])[end]

# Prey
prey = poly!(scene,
       decompose(Point2f0, Circle(Point2f0(0.,0.), preyRadius)),
       color = RGBA(0.9, 0.75, 0.65, 0.5), facealpha = 0.1,
       strokewidth = .25, strokecolor = :black)
preyGut = poly!(scene,
      decompose(Point2f0, Circle(Point2f0(0.,0.), preyRadius-preyMargin)),
      color = RGBA(1., 0.75, 0.75, 0.25), facealpha = 0.1,
      strokewidth = 0.0, strokecolor = RGBA(0.9, 0.75, 0.65, 1.0))
receptorOffColor = [RGBA(0.65, 0.75, 0.40, 1.0) for i in 1:nReceptor]
receptorColor = [RGBA(1.0, 1.0, 0.25, 1.0) for i in 1:nReceptor]

receptor = scatter!(scene, receptorLocation ,
            markersize = 10, color = receptorColor,
            strokecolor = :black, strokewidth = 0.25)[end]

preyStep = [0.0, 0.0]

# predatorField = poly!(scene,
#   lift(s->decompose(Point2f0, Circle(Point2f0(predatorLocation),
#   predatorFieldRadius)), t),
#   color = RGBA(.3, 1.0, .5, .5), strokewidth = 0, strokecolor = :black)




# function moveBacteria(bacteria, dx, dy)
#
#     bacteria

N = 250
#for i in 1:N

record(scene, "test.mp4", 1:N) do i

    global predatorStep = 0.95*predatorStep +
                          0.05*randn(2).*50.0  .- predatorLocation ./ 4000.
    global preyStep = 0.95*preyStep + 0.05*randn(2).*50.
  #  global preyLocation = preyLocation .+ preyStep
    global predatorLocation = predatorLocation .+ predatorStep .+ preyStep
    global bacteriaPos += ones(nBacteria,1)*preyStep'
    bacteria.markersize =
         bacteriaSize.*(sum(bacteriaPos.^2,dims=2) .< matRadius2)[:]
    receptorColor = receptorOffColor[:]

    for j = 1:nReceptor
       range = sqrt( (predatorLocation[1] -receptorLocation[j][1])^2  +
                    (predatorLocation[2] -receptorLocation[j][2])^2 ) -
                     predatorRadius
       if range<0.0 range = 0.0; end
       receptorState[j] = Int(rand()[] < pOpen(range, V))
    end

    receptorColor[findall(x->x==1, receptorState)] .= RGB(1.0, 1.0, 0.0)
    receptor.color[] = receptorColor
    likelihood()
    #LhdPlot.levels[] =  2; #maximum(LikelihoodArray)*[0.1, .5 , .9]
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

# likelihood()
# contour!(x,y,LikelihoodArray,
#           colormap = :Reds,
#           levels = maximum(LikelihoodArray)*[0.1, .5 , .9])
