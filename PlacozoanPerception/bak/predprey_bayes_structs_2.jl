# Predprey_bayes_struct_1
# Placozoan Bayesian particle filter
# with structs

using Makie
using Colors
using OffsetArrays

# Param structure
# physical parameters
struct Param

  # single-cell dipole source
  ρ::Float64    # Resisitivity of seawater Ω.cm
  δ::Float64    # dipole separation cm
  I::Float64    # dipole current A

  # Johnson-Nyquist noise
  kB::Float64       # Bolzmann's constant
  T::Float64        # Kelvin
  Ω::Float64        # receptor impedance Ω
  Δf::Float64       # bandwidth
  σ::Float64        # Johnson-Nyquist noise RMS
end

# Param constructor
function Param()

  ρ = 25.0          # Resisitivity of seawater 25Ω.cm
  δ = 20.e-6*100.   # dipole separation 10μm in cm
  I = 0.1e-12       # dipole current 0.1pA

  # Johnson-Nyquist noise
  kB = 1.38e-23           # Bolzmann's constant
  T = 300.                # Kelvin
  Ω = 20.e6               # receptor impedance
  Δf = 1.0e3              # bandwidth
  σ = sqrt(4.0*kB*T*Δf)   # Johnson-Nyquist noise RMS

  return Param(ρ, δ, I, kB, T, Ω, Δf, σ)
end

# construct parameter structure
param = Param()

# World structure
# contains scene settings and global parameters
struct World
  radius::Int64
  matcolor::RGBA{Float64}
  bgcolor::RGBA{Float64}
  likelihood::OffsetArray     # likelihood given all receptor states
  prior::OffsetArray
  posterior::OffsetArray
end

# World constructor
function World(radius::Int64)
  indx = -radius:radius
  n_indx = length(indx)
  likelihood = OffsetArray(fill(0.0, n_indx,n_indx), indx, indx)
  prior = OffsetArray(fill(0.0, n_indx,n_indx), indx, indx)
  posterior = OffsetArray(fill(0.0, n_indx,n_indx), indx, indx)

  default_matcolor = RGBA(.1, .40, .1, 1.0)
  default_bgcolor  = RGBA(0.0, 0.0, 0.0, 1.0)
  return World(radius, default_matcolor, default_bgcolor,
               likelihood, prior, posterior)
end

# electroreceptor array
# including Bayesian receptive fields (likelihoods for prey proximity)
struct Ereceptor
  N::Int64
  size::Float64  # symbol size for drawing receptors
  x::Array{Float64,1}  # receptor x-coords relative to centre of placozoan
  y::Array{Float64,1}  # receptor y-coords
  Open::Array{Float64,1} # 0/1 for receptor in closed/open state
  openColor::RGB
  closedColor::RGB
  # pOpen[i] is an array containing pre-computed probability pOpen[i][j,k]
  #   that the ith receptor will be in the OPEN state if the nearest edge of
  #   the predator is at the [j,k]th grid point. This is the Bayesian
  #   receptive field of the receptor a.k.a. the likelihood function for
  #   predator proximity given that the receptor is activated (in open state).
  pOpen::Array{OffsetArray,1}
end

# Placozoan structure
struct Placozoan
  radius::Float64
  marginwidth::Float64
  gutradius::Float64
  celldiam::Float64
  x::Array{Float64,1}  # x-ccord of centre
  y::Array{Float64,1}   # y-coord of centre
  # field[i] is pre-computed bio-electric field strength
  #   at distance i μm from edge of placozoan
  field::Array{Float64,1}
  potential::Array{Float64,1}  # in μV
  fieldrange::Int64   # number of elements in field (= max range in μm)
  receptor::Ereceptor
  speed::Array{Float64,1}
  color::RGBA{Float64}
  gutcolor::RGBA{Float64}
  edgecolor::RGB{Float64}
end

# placozoan constructor
# specify size and margin width
# other parameters take default values; located at origin
function Placozoan(radius, margin)
    fieldrange = Int(round(radius*3))
    return Placozoan(radius, margin, radius-margin, 12.0, [0.0], [0.0],
            fill(0.0, fieldrange), fill(0.0, fieldrange), fieldrange,
            [0.0], [],
            RGBA(0.9, 0.75, 0.65, 0.5), RGBA(1., 0.75, 0.75, 0.25),
            RGB(0.0, 0.0, 0.0))
end



# function computes receptor channel Open probability
# as a function of electric field strength
# calibrated to 10% thermal noise-driven open probability for target at infinity
v0 = -param.σ*log(0.1/(1.0-0.1))
pOpenGivenFieldstrength(e) =  1.0./(1 .+ exp.(-(e.-v0)/param.σ))

# function computes single-cell dipole field strength at distance r, in μV/cm
dipoleFieldstrength(r::Float64) = 1.0e6*2π*param.ρ*param.I*param.δ./r.^3

# function computes field strength at distance d μm from edge of body
# due to all dipoles in a placozoan, inserts in placozoan.field
function placozoanFieldstrength(p::Placozoan)
  for a in p.celldiam:p.celldiam:(p.gutradius - p.celldiam)
    n = round(2π*a/p.celldiam)    # number of dipoles in layer
    x = [ a*cos(2π*i/n) for i in 1:n]     # location of dipole
    y = [ a*sin(2π*i/n) for i in 1:n]
    for d in 1:p.fieldrange
      r = sqrt.(((d.+p.radius.-x).^2 + y.^2)).*1.0e-4
      p.field[d] = p.field[d] + sum(dipoleFieldstrength.(r))
    end
    #scatter!(x, y, markersize = cellDiam/16)
  end
end

# compute microvolts across receptor from electric field
function microvoltsFromField(p::Placozoan)
  V = cumsum(p.field)*1.0e-4
  p.potential[:] = V[end].-V
end

function pOpen(d, V)
   i = Int(round(d)) + 1
   if i > length(V)
     i = length(V)
   end
   return pOpenGivenFieldstrength(V[i]*1.0e-6)
 end

# Ereceptor constructor
# creates N receptors centred at (0,0)
# precomputes Bayesian RFs
function Ereceptor(w::World, self_radius, other::Placozoan,
                   N::Int64, displaysize::Float64)
   if floor(N/4)!=N/4
     error("Number of receptors must be a multiple of 4")
   end

   # N receptors equally spaced around edge of placozoan
   x = [self_radius.*(cos(2π*i/N)) for i in 1:N]
   y = [self_radius.*(sin(2π*i/N)) for i in 1:N]
   Open = zeros(N) # initialize receptors closed  (init. state doesn't matter)

   openColor   = RGB(1.0, 1.0, 0.25)
   closedColor = RGB(0.35, 0.45, 0.35)

   # 1d vector containing N offset arrays; ith will contain RF for ith receptor
   Lhd = Array{OffsetArray,1}(undef,N)
   indx = -w.radius:w.radius   # indices for offset array
   n_indx = length(indx)
   for i in 1:N  # for each receptor
     # create an offset array to hold the precomputed RF for that receptor
     Lhd[i] = OffsetArray(fill(0.0, n_indx, n_indx), indx, indx)
     # precompute likelihood (open state probability) for this receptor
    for j in indx
       for k in indx
         Lhd[i][j,k] = pOpen(sqrt((x[i]-j)^2 + (y[i]-k)^2), other.potential)
       end
     end
   end

   return Ereceptor(N, displaysize, x, y,
                    zeros(N), openColor, closedColor, Lhd, self)

end


# precomputes array of likelihoods for nearest edge of predator at
#  grid location (i,j), for each receptor,
#  i.e. Bayesian receptive fields for predator detection.



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


# work starts here

# create world
W = World(500)

# create prey
prey = Placozoan(100.0, 25.0)

# create predator
predator = Placozoan(120.0, 25.0)
predator.speed[] = 1.0
θ = π*rand()[] # Random initial heading (from above)
predator.x[] = (W.radius + predator.radius/2)*cos(θ)
predator.y[] = (W.radius + predator.radius/2)*sin(θ)
Δ = 40.

# create receptors
receptor = Ereceptor(W,prey,predator,4,12.0)

# Receptor parameters
# NB open state probability is computed out to distance maxRange
#    at a finite set of sample points. This is used to pre-compute
#    likelihoods at each mat grid point, for each receptor
# nReceptor = 12  # multiple of 4
# receptorSize = 12
# receptorState = fill(0, nReceptor) # receptorState[i] == 1 if receptor i is active
# maxRange = 3.0*preyRadius  # max sensor range (from centre)
# nRange = Int(maxRange-preyRadius)  # number of sample points in sensor range
# receptorLocation = [preyRadius.*(cos(2π*i/nReceptor), sin(2π*i/nReceptor))
#     for i in 1:nReceptor]  # place receptors around edge of prey
# for j in 1:nReceptor  # initial random receptor states
#   receptorState[j] = Int(rand()[] < 0.1 )
# end

# # Array for pre-computed likelihoods (each mat grid point, each receptor)
# LikelihoodLookup = fill(1.0e-10, nReceptor, Ngrid, Ngrid)
# # Array for likelihood given receptor state
# LikelihoodArray = fill(0.0, Ngrid, Ngrid)

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

# open state probability for source edge at (x,y) receptor at (x0,y0)
function pOpen(iReceptor, x0, y0)
  for i in 1:length(x)
    for j in 1:length(y)
       LikelihoodLookup[iReceptor, i,j] = pOpen(sqrt((x[i]-x0)^2 + (y[j]-y0)^2), V)
     end
   end
 end

# compute likelihood function in world given receptor state
 function likelihood(w::World, receptor::Ereceptor)

    w.likelihood .= 1.0
    for i = 1:receptor.N
      if receptor.Open[i]==1
        w.likelihood .*= receptor.pOpen[i]
      else
        w.likelihood .*= (1.0 .- receptor.pOpen[i])
      end
    end

   # xx = findall(abs.(x).<=preyRadius)
   # yy = findall(abs.(y).<=preyRadius)

    #
    # for i in xx
    #   for j in yy
    #     if x[i]^2 + y[j]^2 < preyRadius^2
    #         LikelihoodArray[i, j] = 0.0
    #       end
    #   end
    # end

    w.likelihood ./= maximum(w.likelihood)

  end

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


# construct sensory particles in prey margin
# by reflectObservationg likelihood sample points through skin
function reflectObservation(likelihoodParticle_xy::Array{Float64,2})
  R = sqrt.(likelihoodParticle_xy[:,1].^2 + likelihoodParticle_xy[:,2].^2)
  r = (preyRadius .- preyMargin*(R.-preyRadius)./(matRadius-preyRadius))::Array{Float64,1}
  #return (r.*xLhdSample./R, r.*yLhdSample./R)
  observationPlot[1] = r.*likelihoodParticle_xy[:,1]./R            # update reflected sample plot
  observationPlot[2] = r.*likelihoodParticle_xy[:,2]./R
end

function reflectBelief(beliefParticle_xy)
  R = sqrt.(beliefParticle_xy[:,1].^2 + beliefParticle_xy[:,2].^2)
  r = (preyRadius .- preyMargin*(R.-preyRadius)./(matRadius-preyRadius))::Array{Float64,1}
  #return (r.*xLhdSample./R, r.*yLhdSample./R)
  beliefPlot[1] = r.*beliefParticle_xy[:,1]./R            # update reflected sample plot
  beliefPlot[2] = r.*beliefParticle_xy[:,2]./R
end

# precompute field strength and voltage from edge of predator
placozoanFieldstrength(predator)   # 1d array
microvoltsFromField(predator)

# time observable
# used to force scene update (nothing depends explicitly on time)
t = Node(0.0)

# construct scene
WorldSize = 2*W.radius+1
scene = Scene(resolution = (WorldSize, WorldSize),
              limits = FRect(-W.radius, -W.radius,WorldSize,WorldSize ),
              show_axis=false, backgroundcolor=W.bgcolor)

# mat is a dark green disc
mat = poly!(scene,
       decompose(Point2f0, Circle(Point2f0(0.,0.), W.radius)),
       color = W.matcolor, strokewidth = 0, strokecolor = :black)

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
predatorShow = poly!(scene,
      lift(s->decompose(Point2f0, Circle(Point2f0(predator.x[], predator.y[]),
      predator.radius)), t),
      color = predator.color, strokewidth = .25, strokecolor = :black)


#
# # for each receptor construct lookup table
# # for likelihood of predator at (x,y) given receptor state.
# # Computed for receptors in first quadrant only,
# # likelihoods for other quadrants are obtained by rotating this 90deg x 3
# Nq = nReceptor ÷ 4  # number of receptors per quadrant
# for iReceptor in 1:Nq
#   pOpen(iReceptor, receptorLocation[iReceptor][1], receptorLocation[iReceptor][2])
# end
# # 2nd-4th quadrants
# for iReceptor in 1:Nq
#   # 2nd
#   for i in 1:Ngrid
#     for j in 1:Ngrid
#       LikelihoodLookup[Nq+iReceptor, i,j] =
#                             LikelihoodLookup[iReceptor, j,Ngrid+1-i]
#       LikelihoodLookup[2*Nq+iReceptor, i,j] =
#                             LikelihoodLookup[iReceptor, Ngrid+1-i,Ngrid+1-j]
#       LikelihoodLookup[3*Nq+iReceptor, i,j] =
#                             LikelihoodLookup[iReceptor, Ngrid + 1 - j,i]
#     end
#   end
#
# end


# likelihood()  # compute likelihood given initial receptor states

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

# duplicate belief particles that collide with observatio n particles
function collision(PParticle, Lparticle)
  #nCollide = 0
  for i in 1:nPosterior_particles
    ix = Int(round(PParticle[i,1] + matRadius))  # x-grid coord ith particle
    for j in 1:nLparticles
      if ix==Lparticle[j,1]  # found matching x-coord
        if Int(round(PParticle[i,2] + matRadius))==Lparticle[j,2] #&y-coord
          ireplace = rand(1:nPosterior_particles)[]  # pick particle to replace
          PParticle[ireplace,:] = PParticle[i,:] + 5.0*randn(2)
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
preyShow = poly!(scene,
       decompose(Point2f0, Circle(Point2f0(0.,0.), prey.radius)),
       color = prey.color, strokewidth = .25, strokecolor = :black)
preyGutShow = poly!(scene,
      decompose(Point2f0, Circle(Point2f0(0.,0.), prey.gutradius)),
      color = prey.gutcolor, strokewidth = 0.0)


receptorShow = scatter!(scene, receptor.x, receptor.y ,
            markersize = receptor.size, color = receptor.openColor,
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

record(scene, "test.mp4", framerate = 24, 1:Nframes) do i

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
          0.1*v*predatorSpeed.*(predatorLocation) ./ d


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
