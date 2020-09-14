# PlacozoanPredatorPrey_01

#import BayesianPlacozoan

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
W = World(500, 1000, 1000)

# create prey
prey = Placozoan(100.0, 25.0)

# create predator
predator = Placozoan(120.0, 25.0)
predator.speed[] = 1.0
θ = π*rand()[] # Random initial heading (from above)
predator.x[] = (W.radius + predator.radius/2)*cos(θ)
predator.y[] = (W.radius + predator.radius/2)*sin(θ)
Δ = 40.


# compute field and potential as a fcn of distance from edge of predator
placozoanFieldstrength(predator)

# compute Bayesian receptive fields for each of prey's receptors
precomputeBayesianRF(W, prey, predator)

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
clock =text!(scene,"t = 0.0s",textsize = 24, color = :white,
     position = (- 0.9*W.radius , -0.9*W.radius))[end]

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




likelihood(W, prey)  # initialize likelihood given initial receptor states
sample_likelihood(W) # sample from normalized likelihood

initialize_posterior()


# plot likelihood particles (samples from likelihood)
LparticlePlot = scatter!(W.Lparticle[:,1], W.Lparticle[:,2],
          color = W.likelycolor, markersize = W.likelysize[],
          strokewidth = 0.1)[end]

# plot projection of likelihood particles into prey margin
# nb this is a dummy plot
# the correct particle locations are inserted before first plot
observationPlot = scatter!(scene,zeros(nLparticles),zero(nLparticles),
      color = :yellow, strokewidth = 0, markersize=2)[end]

# plot projection of posterior particles into prey margin
# nb this is a dummy plot
# the correct particle locations are inserted before first plot
beliefPlot = scatter!(scene,
              zeros(nPosterior_particles), zeros(nPosterior_particles),
              color = :cyan, strokewidth = 0, markersize=2)[end]

#
# #reflectObservation(xParticle, yParticle)
#
# initialize_posterior()
# # plot likelihood particles (sample points)
# PParticlePlot = scatter!(PParticle[:,1], PParticle[:,2],
#           color = postColor, markersize = postSize, strokewidth = 0.1)[end]
#
#
# Prey
preyShow = poly!(scene,
       decompose(Point2f0, Circle(Point2f0(0.,0.), prey.radius)),
       color = prey.color, strokewidth = .25, strokecolor = :black)
preyGutShow = poly!(scene,
      decompose(Point2f0, Circle(Point2f0(0.,0.), prey.gutradius)),
      color = prey.gutcolor, strokewidth = 0.0)


receptorShow = scatter!(scene, prey.receptor.x, prey.receptor.y ,
            markersize = prey.receptor.size, color = prey.receptor.openColor,
            strokecolor = :black, strokewidth = 0.25)[end]

# preyStep = [0.0, 0.0]
# predatorStep = [0.0, 0.0]
# posteriorStep = fill(0.0, nPosterior_particles,2)
#
#
#
# dψ = π/Nframes
# C = [cos(dψ) sin(dψ); -sin(dψ) cos(dψ)]
#
# record(scene, "test.mp4", framerate = 24, 1:Nframes) do i
#
#   global Δ
#   global Nframes
#
#   if i > 0.75*Nframes
#     Δ = Δ + 1000.0/Nframes
#   end
#
#     # predator movement
#     d = sqrt(sum(predatorLocation.^2))  # distance from origin
#     v = sign(preyRadius + predatorRadius + Δ - d)#(distance between edges)-Δ.
#     # pink noise motion in mat frame
#     global predatorStep = 0.8*predatorStep +
#           0.2*randn(2).*predatorSpeed  .+
#           0.1*v*predatorSpeed.*(predatorLocation) ./ d
#
#
#     # prey motion = pink noise in mat frame
#     global preyStep = 0.9*preyStep + 0.1*randn(2).*preySpeed
#
#     # predator location in prey frame
#     global predatorLocation = C*(predatorLocation .+ predatorStep .+ preyStep)
#
#     # # bacteria location in prey frame
#     # global bacteriaPos += ones(nBacteria,1)*preyStep'
#     #
#     # # update bacteria plot
#     # bacteria.markersize =
#     #      bacteriaSize.*(sum(bacteriaPos.^2,dims=2) .< matRadius2)[:]
#
#     updateReceptorState(receptorState)
#
#     likelihood()            # compute likelihood given receptor states
#     sample_likelihood()     # draw random sample from likelihood (indices)
#     # update sample plot
#     likelihoodParticle_xy = -matRadius.+ Lparticle.*sceneWidth/Ngrid
#   #  yy = -matRadius.+ Lparticle[:,2].*sceneWidth/Ngrid
#     LparticlePlot[1] = likelihoodParticle_xy[:,1]
#     LparticlePlot[2] = likelihoodParticle_xy[:,2]
#
#     # xParticle = x[Int.(Lparticle[:,1])]  # convert grid indices to x-y coords
#     # yParticle = y[Int.(Lparticle[:,2])]
#
# #    s = reflectObservation(xParticle, yParticle)  # reflect samples into margin
#     reflectObservation(likelihoodParticle_xy)  # reflect samples into margin
#     # observationPlot[1] = s[1]            # update reflected sample plot
#     # observationPlot[2] = s[2]
#
#     # posterior
#     # pink noise walk (particles mimic predator dynamics)
#     global posteriorStep = 0.95*posteriorStep +
#           0.5*randn(nPosterior_particles,2).*predatorSpeed
#     global PParticle += posteriorStep
#
#    diffusionBarriers()
#
#    # check for collisions between belief particles and observation particles
#    PParticle = collision(PParticle, Lparticle)
#
#     PParticlePlot[1] = PParticle[:,1]
#     PParticlePlot[2] = PParticle[:,2]
#
#     reflectBelief(PParticle)
#
#
#     # level curves of likelihood
#     #LhdPlot.levels[] = maximum(LikelihoodArray)*[0.1, .5 , .9]
#
#     # clock display
#     clock[1] = "t = " * string(floor(t[])) * "s"
#
#     # Node update causes redraw
#     t[] = dt*(i+1)
#
# end


display(scene)
