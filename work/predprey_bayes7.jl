# Bayes 2
# Compute Posterior samples
# 3 Sept 2020

using Makie
using Colors






sceneWidth  = 500.0
matRadius = sceneWidth / 2.0
matRadius2  = matRadius^2
sceneLimits = FRect(-matRadius, -matRadius, sceneWidth, sceneWidth)
#sceneResolution = 1001  # nb must be odd, to ensure grid point @ centre
Ngrid = 501   # grid points in each direction, odd so centre is a grid point
x = LinRange(-matRadius, matRadius, Ngrid)
y = LinRange(-matRadius, matRadius, Ngrid)
dt = 2.5
Nframes = 480

# Prey parameters
preyRadius = 100.
preyMargin = 25.
preySpeed = 0.

# Predator parameters
predatorRadius = 120.   # body radius
predatorMargin = 25.
predatorSpeed = 2.
cellDiam = 12.   # cell radius (used to compute number of dipoles in predator)
θ = π*rand()[] # Random approach heading
predatorLocation = (matRadius+predatorRadius).*(cos(θ),sin(θ)) # initial location
#predatorLocation = (100.0+preyRadius+predatorRadius).*(cos(θ),sin(θ))
Δ = 20.

# Receptor parameters
# NB open state probability is computed out to distance maxRange
#    at a finite set of sample points. This is used to pre-compute
#    likelihoods at each mat grid point, for each receptor
nReceptor = 12  # multiple of 4
receptorSize = 12
receptorState = fill(0, nReceptor) # receptorState[i] == 1 if receptor i is active
maxRange = 3.0*preyRadius  # max sensor range (from centre)
nRange = Int(maxRange-preyRadius)  # number of sample points in sensor range
receptorLocation = [preyRadius.*(cos(2π*i/nReceptor), sin(2π*i/nReceptor))
    for i in 1:nReceptor]  # place receptors around edge of prey
for j in 1:nReceptor  # initial random receptor states
  receptorState[j] = Int(rand()[] < 0.1 )
end

# Array for pre-computed likelihoods (each mat grid point, each receptor)
LikelihoodLookup = fill(1.0e-10, nReceptor, Ngrid, Ngrid)
# Array for likelihood given receptor state
LikelihoodArray = fill(0.0, Ngrid, Ngrid)

# Likelihood particles
nLparticles = 1600   # number of particles in likelihood sample
Lparticle = fill(0, (nLparticles,2)) # grid coords of Likelihood particles
likelyColor = RGB(0.85, 0.65, 0.35)
likleySize = 5

# Posterior density parameters
nPosterior_particles = 800
PParticle = fill(0.0, (nPosterior_particles,2)) # posterior particle locations
postColor = RGB(0.99, 0.35, 0.85)
postSize = 5
priorSD = 75.0

# bacteria parameters
# bacteria are just for show - visualise how prey is moving on mat
bacteriaDensity = 5e-5 # bacteria /um^2
bacteriaColor = RGB(0.5, 0.25, 0.0)
bacteriaRadius = sceneWidth
nBacteria = Int(round(bacteriaDensity*π*bacteriaRadius^2))
bacteriaPos = bacteriaRadius*(rand(nBacteria,2) .- 0.5)
bacteriaSize = 8.0*rand(nBacteria)

# Current dipole parameters
# predator has an array of dipoles generating a proton current
ρ = 25.0   # Resisitivity of seawater 25Ω.cm
δ = 20.e-6*100.  # dipole separation 10μm in cm
I = 0.1e-12     # dipole current 0.1pA

# Johnson-Nyquist noise
# assume receptor cell impedance 20MΩ
kB = 1.38e-23
T = 300.
Ω = 20.e6   # 20MΩ
Δf = 1.0e3   # 1KHz
σ = sqrt(4.0*kB*T*Δf)

# receptor channel open probability as a function of electric field strength
# calibrated to 10% open probability for target at infinity
# (when channel opening is caused by thermal noise)
kT = kB*T
v0 = -σ*log(0.1/(1.0-0.1))
p(e) =  1.0./(1 .+ exp.(-(e.-v0)/σ))

# dipole field strength at distance r from source , μV/cm
E_(r) = 1.0e6*2π*ρ*I*δ./r.^3
# array holds precomputed field strength
# as a function of distance from edge of predator
E = fill(0.0, nRange)

# function returns particle distances from origin
# particle_ij is nParticles x 2, grid coords of particle
function d2o(particle_xy)
  d = fill(0.0, size(particle_xy,1))
  for i in 1:size(particle_xy,1)
    xx = -matRadius + particle_xy[i,1].*sceneWidth/Ngrid
    yy = -matRadius + particle_xy[i,2].*sceneWidth/Ngrid
    d[i] = sqrt(xx^2 + yy^2)
  end
  return d
end





# function to compute field strength at distance d from edge of body
function fieldStrength(sourceRadius, sourceMargin, nRange, E)
  for a in cellDiam:cellDiam:(sourceRadius-cellDiam)
    n = round(2π*a/cellDiam)
    x = [ a*cos(2π*i/n) for i in 1:n]
    y = [ a*sin(2π*i/n) for i in 1:n]
    for d in 1:nRange
      r = sqrt.(((d.+sourceRadius.-x).^2 + y.^2)).*1.0e-4
      E[d] = E[d] + sum(E_(r))
    end
    #scatter!(x, y, markersize = cellDiam/16)
  end
  return E
end

function voltageMicrovolts(d, E)
  V = cumsum(E)*1.0e-4
  return V[end].-V[Int.(round.(d))]
end



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

# construct sensory particles in prey margin
# by reflectObservationg likelihood sample points through skin
function reflectObservation(likelihoodParticle_xy)
  R = sqrt.(likelihoodParticle_xy[:,1].^2 + likelihoodParticle_xy[:,2].^2)
  r = preyRadius .- preyMargin*(R.-preyRadius)./(matRadius-preyRadius)
  #return (r.*xLhdSample./R, r.*yLhdSample./R)
  observationPlot[1] = r.*likelihoodParticle_xy[:,1]./R            # update reflected sample plot
  observationPlot[2] = r.*likelihoodParticle_xy[:,2]./R
end

function reflectBelief(beliefParticle_xy)
  R = sqrt.(beliefParticle_xy[:,1].^2 + beliefParticle_xy[:,2].^2)
  r = preyRadius .- preyMargin*(R.-preyRadius)./(matRadius-preyRadius)
  #return (r.*xLhdSample./R, r.*yLhdSample./R)
  beliefPlot[1] = r.*beliefParticle_xy[:,1]./R            # update reflected sample plot
  beliefPlot[2] = r.*beliefParticle_xy[:,2]./R
end

# precompute field strength and voltage from edge of predator
E = fieldStrength(predatorRadius, predatorMargin, nRange, E)
V = voltageMicrovolts(1:nRange, E)

# time observable
# used to force scene update (nothing depends explicitly on time)
t = Node(0.0)

# construct scene
scene = Scene(resolution = (1000, 1000),
              limits = sceneLimits, show_axis=false,
              backgroundcolor=:black)

# mat is a dark green disc
mat = poly!(scene,
       decompose(Point2f0, Circle(Point2f0(0.,0.), sqrt(matRadius2))),
       color = RGBA(.1, .40, .1, 1.0), strokewidth = 0, strokecolor = :black)

# display nominal time on background
clock =text!(scene,"t = 0.0s",textsize = 12, color = :white,
     position = (- 0.9*matRadius , -0.9*matRadius))[end]

 # scatter bacteria over the mat
 # initial bacteria coords extend off the edge of the mat
 # so that bacteria in fixed locations enter and leave the scene
 # as the prey moves (scene moves in prey frame)
# bacteria = scatter!(scene,
#              lift(s->(bacteriaPos[:,1], bacteriaPos[:,2]), t),
#              color = bacteriaColor, markersize = bacteriaSize,
#              strokewidth = 0, strokecolor = :green)[end]

# predator drawn using lift(..., node)
# (predatorLocation does not depend explicitly on t, but this causes
#  the plot to be updated when the node t changes)
predator = poly!(scene,
      lift(s->decompose(Point2f0, Circle(Point2f0(predatorLocation),
      predatorRadius)), t),
      color = RGBA(.6, 0.6, .5, .75), strokewidth = .25, strokecolor = :black)



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
  for i in 1:Ngrid
    for j in 1:Ngrid
      LikelihoodLookup[Nq+iReceptor, i,j] =
                            LikelihoodLookup[iReceptor, j,Ngrid+1-i]
      LikelihoodLookup[2*Nq+iReceptor, i,j] =
                            LikelihoodLookup[iReceptor, Ngrid+1-i,Ngrid+1-j]
      LikelihoodLookup[3*Nq+iReceptor, i,j] =
                            LikelihoodLookup[iReceptor, Ngrid + 1 - j,i]
    end
  end

end


likelihood()  # compute likelihood given initial receptor states

# Function to sample from normalized likelihood by rejection
function sample_likelihood()

    n = 0
    while n<size(Lparticle,1)
      candidate = rand(1:Ngrid,2)
      if x[candidate[1]]^2 + y[candidate[2]]^2 < matRadius2
        if rand()[]<LikelihoodArray[candidate...]
          n = n + 1
          Lparticle[n, :] = candidate[:]
        end
      end
    end
end

# duplicate belief particles that collide with observation particles
function collision(PParticle, Lparticle)
  #nCollide = 0
  for i in 1:nPosterior_particles
    ix = Int(round(PParticle[i,1] + matRadius))  # x-grid coord ith particle
    for j in 1:nLparticles
      if ix==Lparticle[j,1]  # found matching x-coord
        if Int(round(PParticle[i,2] + matRadius))==Lparticle[j,2] #&y-coord
          ireplace = rand(1:nPosterior_particles)[]  # pick particle to replace
          PParticle[ireplace,:] = PParticle[i,:]
          #nCollide +=1
        end
      end
    end
  end
  #println(nCollide)
  return PParticle
end

sample_likelihood() # sample from normalized likelihood

xParticle = x[Int.(Lparticle[:,1])]  # convert grid indices to x-y coords
yParticle = y[Int.(Lparticle[:,2])]
# plot likelihood particles (sample points)
LparticlePlot = scatter!(xParticle, yParticle,
          color = likelyColor, markersize = likleySize, strokewidth = 0.1)[end]

# project likelihood particles into prey margin and plot

observationPlot = scatter!(scene,xParticle, yParticle,
      color = :yellow, strokewidth = 0, markersize=2)[end]

beliefPlot = scatter!(scene,
              rand(nPosterior_particles), rand(nPosterior_particles),
            color = :cyan, strokewidth = 0, markersize=2)[end]


#reflectObservation(xParticle, yParticle)
# initialize posterior samples
# truncated Gaussian distribution of distance from mat edge
# (i.e. diffusion from edge with absorbing barrier at prey)
function initialize_posterior()

  nP = 0
  while nP < nPosterior_particles
    ϕ = 2.0*π*rand(1)[]
    β = 1.0e12
    while β > (matRadius-preyRadius)
      β = priorSD*abs(randn(1)[])
    end
    #candidate = -matRadius .+ sceneWidth.*rand(2)  # random point in scene
    # d = sqrt(candidate[1]^2 + candidate[2]^2) # candidate distance from origin
    # if (d>preyRadius) & (d<matRadius)
      nP = nP+1
      PParticle[nP,:] =  (matRadius-β).*[cos(ϕ), sin(ϕ)]
    # end
  end
end

initialize_posterior()
# plot likelihood particles (sample points)
PParticlePlot = scatter!(PParticle[:,1], PParticle[:,2],
          color = postColor, markersize = postSize, strokewidth = 0.1)[end]


# Prey
prey = poly!(scene,
       decompose(Point2f0, Circle(Point2f0(0.,0.), preyRadius)),
       color = RGBA(0.9, 0.75, 0.65, 0.5), facealpha = 0.1,
       strokewidth = .25, strokecolor = :black)
preyGut = poly!(scene,
      decompose(Point2f0, Circle(Point2f0(0.,0.), preyRadius-preyMargin)),
      color = RGBA(1., 0.75, 0.75, 0.25), facealpha = 0.1,
      strokewidth = 0.0, strokecolor = RGBA(0.9, 0.75, 0.65, 1.0))
receptorOffColor = [RGBA(0.35, 0.45, 0.35, 1.0) for i in 1:nReceptor]
receptorColor = [RGBA(1.0, 1.0, 0.25, 1.0) for i in 1:nReceptor]

receptor = scatter!(scene, receptorLocation ,
            markersize = receptorSize, color = receptorColor,
            strokecolor = :black, strokewidth = 0.25)[end]

preyStep = [0.0, 0.0]
predatorStep = [0.0, 0.0]
posteriorStep = fill(0.0, nPosterior_particles,2)


function updateReceptorState(receptorState)

  # turn all receptors off
  receptorColor = receptorOffColor[:]

  # calculate receptor states
  for j = 1:length(receptorState)
     range = sqrt( (predatorLocation[1] -receptorLocation[j][1])^2  +
                  (predatorLocation[2] -receptorLocation[j][2])^2 ) -
                   predatorRadius
     if range<0.0 range = 0.0; end
     receptorState[j] = Int(rand()[] < pOpen(range, V))
  end

  # update receptor state display
  receptorColor[findall(x->x==1, receptorState)] .= RGB(1.0, 1.0, 0.0)
  receptor.color[] = receptorColor
end

function diffusionBarriers()
  # posterior particle diffusion barriers at edge of prey and of mat
  for j in 1:nPosterior_particles
    d = sqrt(PParticle[j,1]^2+PParticle[j,2]^2)
      if d>matRadius # edge of mat is reflectObservationg barrier
         PParticle[j,:] = matRadius.*PParticle[j,:]./d
         posteriorStep[j,:] = [0.0, 0.0]
      end
      if d<preyRadius  # edge of prey is absorbing barrier
         ϕ = 2.0*π*rand(1)[]
         PParticle[j,:] = matRadius.*[cos(ϕ), sin(ϕ)]
         posteriorStep[j,:] = [0.0, 0.0]
      end
  end
end

dψ = π/Nframes
C = [cos(dψ) sin(dψ); -sin(dψ) cos(dψ)]

record(scene, "test.mp4", framerate = 12, 1:Nframes) do i

  global Δ
  global Nframes

  if i > 0.75*Nframes
    Δ = Δ + 1000.0/Nframes
  end

    # predator movement
    d = sqrt(sum(predatorLocation.^2))  # distance from origin
    v = sign(preyRadius + predatorRadius + Δ - d)#(distance between edges)-Δ.
    # pink noise motion in mat frame
    global predatorStep = 0.8*predatorStep +
          0.2*randn(2).*predatorSpeed  .+
          0.05*v*predatorSpeed.*(predatorLocation) ./ d


    # prey motion = pink noise in mat frame
    global preyStep = 0.9*preyStep + 0.1*randn(2).*preySpeed

    # predator location in prey frame
    global predatorLocation = C*(predatorLocation .+ predatorStep .+ preyStep)

    # # bacteria location in prey frame
    # global bacteriaPos += ones(nBacteria,1)*preyStep'
    #
    # # update bacteria plot
    # bacteria.markersize =
    #      bacteriaSize.*(sum(bacteriaPos.^2,dims=2) .< matRadius2)[:]

    updateReceptorState(receptorState)

    likelihood()            # compute likelihood given receptor states
    sample_likelihood()     # draw random sample from likelihood (indices)
    # update sample plot
    likelihoodParticle_xy = -matRadius.+ Lparticle.*sceneWidth/Ngrid
  #  yy = -matRadius.+ Lparticle[:,2].*sceneWidth/Ngrid
    LparticlePlot[1] = likelihoodParticle_xy[:,1]
    LparticlePlot[2] = likelihoodParticle_xy[:,2]

    # xParticle = x[Int.(Lparticle[:,1])]  # convert grid indices to x-y coords
    # yParticle = y[Int.(Lparticle[:,2])]

#    s = reflectObservation(xParticle, yParticle)  # reflect samples into margin
    reflectObservation(likelihoodParticle_xy)  # reflect samples into margin
    # observationPlot[1] = s[1]            # update reflected sample plot
    # observationPlot[2] = s[2]

    # posterior
    # pink noise walk (particles mimic predator dynamics)
    global posteriorStep = 0.95*posteriorStep +
          0.5*randn(nPosterior_particles,2).*predatorSpeed
    global PParticle += posteriorStep

   diffusionBarriers()

   # check for collisions between belief particles and observation particles
   PParticle = collision(PParticle, Lparticle)

    PParticlePlot[1] = PParticle[:,1]
    PParticlePlot[2] = PParticle[:,2]

    reflectBelief(PParticle)


    # level curves of likelihood
    #LhdPlot.levels[] = maximum(LikelihoodArray)*[0.1, .5 , .9]

    # clock display
    clock[1] = "t = " * string(floor(t[])) * "s"

    # Node update causes redraw
    t[] = dt*(i+1)

end


display(scene)
