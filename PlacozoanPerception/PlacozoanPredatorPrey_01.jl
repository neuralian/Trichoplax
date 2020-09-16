# PlacozoanPredatorPrey_01

#import BayesianPlacozoan

# sceneWidth  = 500.0
# matRadius = sceneWidth / 2.0
# matRadius2  = matRadius^2
# sceneLimits = FRect(-matRadius, -matRadius, sceneWidth, sceneWidth)
# #sceneResolution = 1001  # nb must be odd, to ensure grid point @ centre
# Ngrid = 501   # grid points in each direction, odd so centre is a grid point
# x = LinRange(-matRadius, matRadius, Ngrid)
# y = LinRange(-matRadius, matRadius, Ngrid)

# Simulation parameters
nFrames = 500
dt = 1.0

# World/physical parameters
mat_radius = 400
n_likelihood_particles = 500
n_posterior_particles = 500

# placozoan parameters
prey_radius = 80.0
prey_margin = 25.0
Nreceptors = 12

predator_radius = 90.0
predator_margin = 10.0
predator_speed = 1.0

approach_Δ = 40.0   # proximity at which predator stops approaching prey

# create world
W = World(nFrames, mat_radius,
           n_likelihood_particles, n_posterior_particles, approach_Δ)

# create prey
prey = Placozoan(prey_radius, prey_margin, Nreceptors)

# create predator
predator = Placozoan(predator_radius, predator_margin, 0,
                     RGBA(.25, 0.1, 0.1, 1.0),
                     RGBA(.45, 0.1, 0.1, 0.25),
                     RGB(.95, 0.1, 0.1) )
predator.speed[] = predator_speed
θ = π*rand()[] # Random initial heading (from above)
predator.x[] = (mat_radius + predator_radius/2)*cos(θ)
predator.y[] = (mat_radius + predator_radius/2)*sin(θ)
# predator.edgecolor[:] = RGB(.45, 0.1, 0.1)
# predator.color[:] = RGBA(.45, 0.1, 0.1, 0.25)

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

# # Likelihood particles
# nLparticles = 1600   # number of particles in likelihood sample
# Lparticle = fill(0, (nLparticles,2)) # grid coords of Likelihood particles
# likelyColor = RGB(0.85, 0.65, 0.35)
# likleySize = 5
#
# # Posterior density parameters
# nPosterior_particles = 800
# PParticle = fill(0.0, (nPosterior_particles,2)) # posterior particle locations
# postColor = RGB(0.99, 0.35, 0.85)
# postSize = 5
# priorSD = 75.0

# # bacteria parameters
# # bacteria are just for show - visualise how prey is moving on mat
# bacteriaDensity = 5e-5 # bacteria /um^2
# bacteriaColor = RGB(0.5, 0.25, 0.0)
# bacteriaRadius = sceneWidth
# nBacteria = Int(round(bacteriaDensity*π*bacteriaRadius^2))
# bacteriaPos = bacteriaRadius*(rand(nBacteria,2) .- 0.5)
# bacteriaSize = 8.0*rand(nBacteria)

# time observable
# used to force scene update (nothing depends explicitly on time)
t = Node(0.0)

# construct scene
WorldSize = 2*mat_radius+1
scene = Scene(resolution = (WorldSize, WorldSize),
              limits = FRect(-W.radius, -W.radius,WorldSize,WorldSize ),
              show_axis=false, backgroundcolor=W.bgcolor)

# mat is a dark green disc
mat_plt = poly!(scene,
       decompose(Point2f0, Circle(Point2f0(0.,0.), W.radius)),
       color = W.matcolor, strokewidth = 0, strokecolor = :black)

# display nominal time on background
clock_plt =text!(scene,"t = 0.0s",textsize = 24, color = :white,
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
predator_plt = poly!(scene,
      lift(s->decompose(Point2f0, Circle(Point2f0(predator.x[], predator.y[]),
      predator.radius)), t),
      color = predator.color, strokecolor = predator.edgecolor,
      strokewidth = .5)

likelihood(W, prey)  # initialize likelihood given initial receptor states
sample_likelihood(W) # sample from normalized likelihood



# plot likelihood particles (samples from likelihood)
Lparticle_plt = scatter!(W.Lparticle[:,1], W.Lparticle[:,2],
          color = W.likelycolor, markersize = W.likelysize[],
          strokewidth = 0.1)[end]

# plot projection of likelihood particles into prey margin
# nb this is a dummy plot
# the correct particle locations are inserted before first plot
observation_plt = scatter!(scene,zeros(W.nLparticles),zeros(W.nLparticles),
      color = :yellow, strokewidth = 0, markersize=2)[end]

initialize_posterior(W,prey)
Pparticle_plt = scatter!(W.Pparticle[:,1], W.Pparticle[:,2],
          color = W.postcolor, markersize = W.postsize[], strokewidth = 0.1)[end]


# plot projection of posterior particles into prey margin
# nb this is a dummy plot
# the correct particle locations are inserted before first plot
belief_plt = scatter!(scene,
              zeros(W.nPparticles), zeros(W.nPparticles),
              color = :cyan, strokewidth = 0, markersize=2)[end]

# Prey
prey_plt = poly!(scene,
       decompose(Point2f0, Circle(Point2f0(0.,0.), prey.radius)),
       color = prey.color, strokewidth = .25, strokecolor = :black)
preyGut_plt = poly!(scene,
      decompose(Point2f0, Circle(Point2f0(0.,0.), prey.gutradius)),
      color = prey.gutcolor, strokewidth = 0.0)


receptor_plt = scatter!(scene, prey.receptor.x, prey.receptor.y ,
            markersize = prey.receptor.size,
            color = [prey.receptor.openColor for i in 1:prey.receptor.N],
            strokecolor = :black, strokewidth = 0.25)[end]


record(scene, "test.mp4", framerate = 24, 1:W.nFrames) do i

    # cause predator to drift away from prey in last 25% of animation
    if i > 0.75*W.nFrames
      W.Δ[] += 1000.0/W.nFrames
    end

    # predator random walk to within Δ of prey
    # small bias velocity added for drift towards prey & clockwise orbit
    stalk(W, predator, prey)

    # prey receptors respond to predator electric field
    updateReceptors(prey, predator)
    # set color of each receptor, indicating open or closed state
    receptorColor = [prey.receptor.closedColor  for j in 1:prey.receptor.N]
    receptorColor[findall(x->x==1, prey.receptor.state)] .=
         prey.receptor.openColor
    receptor_plt.color[] = receptorColor

    # prey sensory observations (particles released by active sensors)
    likelihood(W, prey)      # likelihood given receptor states
    sample_likelihood(W)     # random sample from likelihood
    Lparticle_plt[1] = W.Lparticle[:,1]   # update likelihood particle plot
    Lparticle_plt[2] = W.Lparticle[:,2]

    observation = reflectObservation(W, prey) # reflect samples into margin
    observation_plt[1] = observation[:,1]     # update observation particle plot
    observation_plt[2] = observation[:,2]

    # predict posterior using predator model
    posteriorPredict(W, predator)
    bayesBelief(W)
    diffusionBoundary(W, prey) # stop particles diffusing out of the arena
    Pparticle_plt[1] = W.Pparticle[:,1]  # update posterior particle plot
    Pparticle_plt[2] = W.Pparticle[:,2]



#    # check for collisions between belief particles and observation particles
#    PParticle = collision(PParticle, Lparticle)
#
#     PParticlePlot[1] = PParticle[:,1]
#     PParticlePlot[2] = PParticle[:,2]
#
#     reflectBelief(PParticle)


    # level curves of likelihood
    #LhdPlot.levels[] = maximum(LikelihoodArray)*[0.1, .5 , .9]

    # clock display
    clock_plt[1] = "t = " * string(floor(t[])) * "s"

    # Node update causes redraw
    t[] = dt*(i+1)

end
