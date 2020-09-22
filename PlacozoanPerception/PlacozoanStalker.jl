# PlacozoanPredatorPrey

#import BayesianPlacozoan

# simulation parameters
nFrames = 1200            # number of animation frames
mat_radius = 365
approach_Δ = 100.0         # predator closest approach distance
dt = 1.00

# construct observer
priorSD = 25.0
posteriorSD = 100.0
n_likelihood_particles = 2500
n_prior_particles = 500
n_posterior_particles = 500
# observer = Observer(mat_radius,
#             n_likelihood_particles, n_prior_particles, n_posterior_particles,
#             approach_Δ)

# construct prey
prey_radius = 100
prey_margin = 25
Nreceptors = 32
prey_fieldrange = 0   # no field
prey = Placozoan(prey_radius, prey_margin, prey_fieldrange,
                  Nreceptors, sizeof_receptor, mat_radius,
                  n_likelihood_particles, n_prior_particles,
                  n_posterior_particles, priorSD, posteriorSD)

# construct predator
# nb has dummy observer
predator_radius = 120
predator_margin = 0
predator_speed = 0.4
predator_fieldrange = mat_radius
predator = Placozoan(predator_radius, predator_margin, predator_fieldrange,
                     RGBA(.25, 0.1, 0.1, 1.0),
                     RGBA(.45, 0.1, 0.1, 0.25),
                     RGB(.95, 0.1, 0.1) )
predator.speed[] = predator_speed
θ = π*rand()[] # Random initial heading (from above)
predator.x[] = (mat_radius + predator_radius)*cos(θ)
predator.y[] = (mat_radius + predator_radius)*sin(θ)


# compute field and potential as a fcn of distance from edge of predator
placozoanFieldstrength!(predator)

# compute Bayesian receptive fields for each of prey's receptors
precomputeBayesianRF(prey, predator)

# time observable
# used to force scene update (nothing depends explicitly on time)
t = Node(0.0)

# construct scene
WorldSize = 2*mat_radius+1
scene = Scene(resolution = (WorldSize, WorldSize),
              limits = FRect(-mat_radius, -mat_radius ,WorldSize, WorldSize ),
              show_axis=false, backgroundcolor = colour_background)

# mat is a dark green disc
mat_plt = poly!(scene,
       decompose(Point2f0, Circle(Point2f0(0.,0.), mat_radius)),
       color = colour_mat, strokewidth = 0, strokecolor = :black)

# display nominal time on background
clock_plt =text!(scene,"t = 0.0s",textsize = 24, color = :white,
     position = (- 0.925*mat_radius , -0.95*mat_radius))[end]

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

likelihood(prey)  # initialize likelihood given initial receptor states
sample_likelihood(prey) # sample from normalized likelihood



# plot likelihood particles (samples from likelihood)
Lparticle_plt = scatter!(
          prey.observer.Lparticle[:,1], prey.observer.Lparticle[:,2],
          color =:yellow, markersize = size_likelihood,
          strokewidth = 0.1)[end]

# plot projection of likelihood particles into prey margin
# nb this is a dummy plot
# the correct particle locations are inserted before first plot
observation_plt = scatter!(scene,
    zeros(prey.observer.nLparticles),zeros(prey.observer.nLparticles),
      color = :yellow, strokewidth = 0, markersize=size_observation)[end]

# initialize_prior_Gaussian(prey)
# Pparticle_plt = scatter!(
#           prey.observer.Pparticle[:,1], prey.observer.Pparticle[:,2],
#           color = colour_prior,
#           markersize = size_prior, strokewidth = 0)[end]

initialize_posterior_Gaussian(prey)
Bparticle_plt = scatter!(
          prey.observer.Bparticle[:,1], prey.observer.Bparticle[:,2],
          color = colour_posterior,
          markersize = size_posterior, strokewidth = 0.1)[end]

# plot projection of posterior particles into prey margin
# nb this is a dummy plot
# the correct particle locations are inserted before first plot
belief_plt = scatter!(scene,
            zeros(prey.observer.nPparticles), zeros(prey.observer.nPparticles),
            color = colour_posterior, strokewidth = 0, markersize=size_belief)[end]

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

#for i in 1:10
record(scene, "PlacozoanPerception.mp4", framerate = 24, 1:nFrames) do i

    # #cause predator to drift away from prey in last 25% of animation
    # if i > 0.85*W.nFrames
    #   W.Δ[] += 400.0/W.nFrames
    # end

    # predator random walk to within Δ of prey
    # small bias velocity added for drift towards prey & clockwise orbit
    stalk(predator, prey, approach_Δ)

    # prey receptors respond to predator electric field
    updateReceptors(prey, predator)
    # set color of each receptor, indicating open or closed state
    receptorColor = [prey.receptor.closedColor  for j in 1:prey.receptor.N]
    receptorColor[findall(x->x==1, prey.receptor.state)] .=
         prey.receptor.openColor
    receptor_plt.color[] = receptorColor

    # prey sensory observations (particles released by active sensors)
    likelihood(prey)      # likelihood given receptor states
    sample_likelihood(prey)     # random sample from likelihood




    # predict posterior using predator model
#    posteriorPredict(W, predator)
    #bayesBelief(W)
    bayesUpdate(prey)


    # steadyPrior(prey) # stop particles diffusing out of the arena


    (observation, belief) = reflect(prey) # reflect samples into margin


    Lparticle_plt[1] = prey.observer.Lparticle[:,1]   # update likelihood particle plot
    Lparticle_plt[2] = prey.observer.Lparticle[:,2]

    observation_plt[1] = observation[:,1]     # update observation particle plot
    observation_plt[2] = observation[:,2]

    # Pparticle_plt[1] = prey.observer.Pparticle[:,1]  # update posterior particle plot
    # Pparticle_plt[2] = prey.observer.Pparticle[:,2]

    Bparticle_plt[1] = prey.observer.Bparticle[:,1]  # update posterior particle plot
    Bparticle_plt[2] = prey.observer.Bparticle[:,2]

    belief_plt[1] = belief[:,1]     # update observation particle plot
    belief_plt[2] = belief[:,2]

    # level curves of likelihood
    #LhdPlot.levels[] = maximum(LikelihoodArray)*[0.1, .5 , .9]

    # clock display
    clock_plt[1] = "t = " * string(floor(t[])) * "s"

    # Node update causes redraw
    t[] = dt*(i+1)

end
