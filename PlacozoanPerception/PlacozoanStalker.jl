# PlacozoanStalker

using AbstractPlotting.MakieLayout
using AbstractPlotting

#import BayesianPlacozoan

# choose display output format: "particles", "arrays" or "all"
PLOT_OUTPUT = "particles"
# PLOT_OUTPUT = "arrays"
# PLOT_OUTPUT = "all"

# simulation parameters
nFrames = 100            # number of animation frames
mat_radius = 400
approach_Δ = 25.0         # predator closest approach distance
dt = 1.00

# construct observer
priormean = 300.
priorsd = 25.0
posteriorSD = 100.0
n_likelihood_particles = 5000
#n_prior_particles = 500
n_posterior_particles = 2500
# observer = Observer(mat_radius,
#             n_likelihood_particles, n_prior_particles, n_posterior_particles,
#             approach_Δ)

# construct prey
prey_radius = 200
prey_margin = 50
Nreceptors = 24
prey_fieldrange = 0   # no field
prey = Placozoan(prey_radius, prey_margin, prey_fieldrange,
                  Nreceptors, sizeof_receptor, mat_radius,
                  n_likelihood_particles, n_posterior_particles,
                  priormean, priorsd)

# construct predator
# nb has dummy observer
predator_radius = 225
predator_margin = 0
predator_speed = 0.6
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

likelihood(prey)  # initialize likelihood given initial receptor states
sample_likelihood(prey) # sample from normalized likelihood
initialize_posterior_particles_Gaussian(prey)
initialize_posterior_array_Gaussian(prey)

# time observable
# used to force scene update (nothing depends explicitly on time)
t = Node(0.0)

# if NO_PLOT == false

# construct scene
WorldSize = 2*mat_radius+1
if PLOT_OUTPUT == "particles"
    scene, layout = layoutscene(resolution = (WorldSize, WorldSize))
    left_panel = layout[1,1] = LAxis(scene, backgroundcolor = colour_background)
    # clock_panel = LAxis(scene,
    #      bbox = BBox(100 , 100, 100, 100) )
    hidespines!(left_panel)
    hidedecorations!(left_panel)
    # hidespines!(clock_panel)
    # hidedecorations!(clock_panel)
end

if PLOT_OUTPUT == "arrays"
    scene, panel = layoutscene(resolution = (2*WorldSize + 20, WorldSize))
    panel[1,1] = LAxis(scene)
    panel[1,2] = LAxis(scene)
    colsize!(panel, 1, fixed(WorldSize))
    colsize!(panel, 2, fixed(WorldSize))
              # limits = FRect(-mat_radius, -mat_radius ,WorldSize, WorldSize ),
              # show_axis=false, backgroundcolor = colour_background)
end


# mat is a dark green disc
mat_plt = poly!(left_panel,
       decompose(Point2f0, Circle(Point2f0(0.,0.), mat_radius)),
       color = colour_mat, strokewidth = 0, strokecolor = :black)

# display nominal time on background
clock_plt =LText(scene,"              t = 0.0s", color = :white)

# predator drawn using lift(..., node)
# (predatorLocation does not depend explicitly on t, but this causes
#  the plot to be updated when the node t changes)
predator_plt = poly!(left_panel,
      lift(s->decompose(Point2f0, Circle(Point2f0(predator.x[], predator.y[]),
      predator.radius)), t),
      color = predator.color, strokecolor = predator.edgecolor,
      strokewidth = .5)



# plot likelihood particles (samples from likelihood)
Lparticle_plt = scatter!(left_panel,
          prey.observer.Lparticle[:,1], prey.observer.Lparticle[:,2],
          color =:yellow, markersize = size_likelihood,
          strokewidth = 0.1)

# plot projection of likelihood particles into prey margin
# nb this is a dummy plot
# the correct particle locations are inserted before first plot
observation_plt = scatter!(left_panel,
    zeros(prey.observer.nLparticles),zeros(prey.observer.nLparticles),
      color = :yellow, strokewidth = 0, markersize=size_observation)

Bparticle_plt = scatter!(left_panel,
          prey.observer.Bparticle[:,1], prey.observer.Bparticle[:,2],
          color = colour_posterior,
          markersize = size_posterior, strokewidth = 0.1)

# plot projection of posterior particles into prey margin
# nb this is a dummy plot
# the correct particle locations are inserted before first plot
belief_plt = scatter!(left_panel,
            zeros(prey.observer.nBparticles), zeros(prey.observer.nBparticles),
            color = colour_posterior, strokewidth = 0, markersize=size_belief)

# Prey
prey_plt = poly!(left_panel,
       decompose(Point2f0, Circle(Point2f0(0.,0.), prey.radius)),
       color = prey.color, strokewidth = 1, strokecolor = RGB(.5, .75, .85))
preyGut_plt = poly!(left_panel,
      decompose(Point2f0, Circle(Point2f0(0.,0.), prey.gutradius)),
      color = prey.gutcolor, strokewidth = 0.0)


receptor_plt = scatter!(left_panel, prey.receptor.x, prey.receptor.y ,
            markersize = prey.receptor.size,
            color = [prey.receptor.openColor for i in 1:prey.receptor.N],
            strokecolor = :black, strokewidth = 0.25)

# reset axis limits (have been auto-adjusted by MakieLayout)
xlims!(left_panel, -mat_radius, mat_radius)
ylims!(left_panel, -mat_radius, mat_radius)


record(scene, "PlacozoanPerception.mp4", framerate = 24, 1:nFrames) do i

    # predator random walk to within Δ of prey
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
    bayesParticleUpdate(prey)
    bayesArrayUpdate(prey)

    (observation, belief) = reflect(prey) # reflect samples into margin

    Lparticle_plt[1] = prey.observer.Lparticle[:,1]   # update likelihood particle plot
    Lparticle_plt[2] = prey.observer.Lparticle[:,2]

    observation_plt[1] = observation[:,1]     # update observation particle plot
    observation_plt[2] = observation[:,2]

    Bparticle_plt[1] = prey.observer.Bparticle[:,1]  # update posterior particle plot
    Bparticle_plt[2] = prey.observer.Bparticle[:,2]

    belief_plt[1] = belief[:,1]     # update observation particle plot
    belief_plt[2] = belief[:,2]

    # level curves of likelihood
    #LhdPlot.levels[] = maximum(LikelihoodArray)*[0.1, .5 , .9]

    # clock display
    clock_plt.text = "              t = " * string(floor(t[])) * "s"

    # Node update causes redraw
    t[] = dt*(i+1)

end

# else  # NO_PLOT == true
#     for i in 1:nFrames
#         stalk(predator, prey, approach_Δ)
#         updateReceptors(prey, predator)
#         likelihood(prey)      # likelihood given receptor states
#         sample_likelihood(prey)     # random sample from likelihood
#         bayesParticleUpdate(prey)
#     end
# end
