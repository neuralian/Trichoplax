# PlacozoanStalker
# function main(nFrames=100)

using TickTock
# import BayesianPlacozoan
# function main()
    
# choose display output
PLOT_EXT_PARTICLES = true
PLOT_INT_PARTICLES = true
PLOT_ARRAYS = true

# simulation parameters
nReps = 2
nFrames = 180        # number of animation frames
mat_radius = 400
approach_Δ = 16.0         # predator closest approach distance
dt = 1.00

    # arrays to record predator location and observer MAP estimate

# construct observer
priorDensity = 0.002
posteriorSD = 100.0
n_likelihood_particles = 5000
# n_prior_particles = 500
n_posterior_particles = 2500
# observer = Observer(mat_radius,
#             n_likelihood_particles, n_prior_particles, n_posterior_particles,
#             approach_Δ)

# construct prey
prey_radius = 120
prey_margin = 30
Nreceptors = 16
Ncrystals = 8
prey_fieldrange = 0   # no field
prey = Placozoan(prey_radius, prey_margin, prey_fieldrange,
                  Nreceptors, sizeof_receptor, mat_radius,
                  Ncrystals, sizeof_crystal, mat_radius,
                  n_likelihood_particles, n_posterior_particles,
                  priorDensity, nFrames)

# construct predator
# nb has dummy observer
predator_radius = 150
predator_margin = 0
predator_speed = 2.0
predator_fieldrange = mat_radius
predator = Placozoan(predator_radius, predator_margin, predator_fieldrange,
                     RGBA(.25, 0.1, 0.1, 1.0),
                     RGBA(.45, 0.1, 0.1, 0.25),
                     RGB(.95, 0.1, 0.1) )
predator.speed[] = predator_speed
θ = π * rand()[] # Random initial heading (from above)
predator.x[] = (mat_radius + 0.5 * predator_radius) * cos(θ)
predator.y[] = (mat_radius + 0.5 * predator_radius) * sin(θ)


# compute field and potential as a fcn of distance from edge of predator
placozoanFieldstrength!(predator)

# pre-compute Bayesian receptive fields for each of prey's receptors
Ereceptor_RF(prey, predator)
Vreceptor_RF(prey)

# initialize observers
likelihood(prey)           # initialize likelihood given initial receptor states
sample_likelihood(prey)    # sample from normalized likelihood
initialize_particles(prey) # draw initial sample from prior
initialize_prior(prey)     # initialize numerical Bayesian prior

# time observable
# used to force scene update (nothing depends explicitly on time)
t = Node(0.0)


# construct scene
WorldSize = 2 * mat_radius + 1
mask = zeros(WorldSize, WorldSize)  # mat mask
for i in 1:WorldSize
    for j in 1:WorldSize
        if ((i - mat_radius - 1)^2 + (j - mat_radius - 1)^2) <= mat_radius^2
            mask[i,j] = 1.0
        end
    end
end
if !PLOT_ARRAYS    # not plotting likelihoods or posterior
    scene, layout = layoutscene(resolution=(WorldSize, WorldSize))
    left_panel =
        layout[1, 1] = LAxis(
            scene,
            title="Placozoan",
            backgroundcolor=colour_background,
        )
    # clock_panel = LAxis(scene,
    #      bbox = BBox(100 , 100, 100, 100) )
    hidespines!(left_panel)
    hidedecorations!(left_panel)
else
    scene, layout = layoutscene(resolution=(3 * .5 * WorldSize + 40, .5 * WorldSize + 20))
    shim = layout[1, 1] = LAxis(scene)
    left_panel =
        layout[1, 2] = LAxis(
            scene,
            title="Placozoan",
            backgroundcolor=colour_background,
        )
    middle_panel = layout[1, 3] = LAxis(scene, title="Likelihood")
    right_panel = layout[1, 4] = LAxis(scene, title="Posterior")
    colsize!(layout, 1, Relative(0.04))
    colsize!(layout, 2, Relative(0.32))
    colsize!(layout, 3, Relative(0.32))
    colsize!(layout, 4, Relative(0.32))
    hidespines!(shim)
    hidedecorations!(shim)
    hidespines!(left_panel)
    hidedecorations!(left_panel)
    hidespines!(middle_panel)
    hidedecorations!(middle_panel)
    hidespines!(right_panel)
    hidedecorations!(right_panel)
end


# mat is a dark green disc in left panel (always present)
mat_plt = poly!(left_panel,
       decompose(Point2f0, Circle(Point2f0(0., 0.), mat_radius)),
       color=colour_mat, strokewidth=0, strokecolor=:black)

if PLOT_ARRAYS
    mat_middle_plt = poly!(
        middle_panel,
        decompose(
            Point2f0,
            Circle(Point2f0(mat_radius, mat_radius), mat_radius),
        ),
        color=RGBA(0.0, 0.0, 0.0, 0.0),
        strokewidth=0.5,
        strokecolor=RGB(0.75, 0.75, 0.75),
    )
    mat_right_plt = poly!(
        right_panel,
        decompose(
            Point2f0,
            Circle(Point2f0(mat_radius, mat_radius), mat_radius),
        ),
        color=RGBA(0.0, 0.0, 0.0, 0.0),
        strokewidth=0.5,
        strokecolor=RGB(0.75, 0.75, 0.75),
    )
end


# display nominal time on background
clock_plt = LText(scene,
             "                                                  t = 0.0s",
             color=:white)

# predator drawn using lift(..., node)
# (predatorLocation does not depend explicitly on t, but this causes
#  the plot to be updated when the node t changes)
predator_plt = poly!(left_panel,
      lift(s -> decompose(Point2f0, Circle(Point2f0(predator.x[], predator.y[]),
      predator.radius)), t),
      color=predator.color, strokecolor=predator.edgecolor,
      strokewidth=.5)

if PLOT_EXT_PARTICLES

    # plot likelihood particles (samples from likelihood)
    Lparticle_plt = scatter!(
        left_panel,
        prey.observer.Lparticle[:, 1],
        prey.observer.Lparticle[:, 2],
        color=:yellow,
        markersize=size_likelihood,
        strokewidth=0.1,
    )

    # plot posterior particles
    Bparticle_plt = scatter!(
        left_panel,
        prey.observer.Bparticle[:, 1],
        prey.observer.Bparticle[:, 2],
        color=colour_posterior,
        markersize=size_posterior,
        strokewidth=0.1,
    )

end # plot external particles

if PLOT_INT_PARTICLES

    # plot projection of likelihood particles into prey margin
    # nb this is a dummy plot
    # the correct particle locations are inserted before first plot
    observation_plt = scatter!(
        left_panel,
        zeros(prey.observer.nLparticles),
        zeros(prey.observer.nLparticles),
        color=:yellow,
        strokewidth=0,
        markersize=size_observation,
    )


    # plot projection of posterior particles into prey margin
    # nb this is a dummy plot
    # the correct particle locations are inserted before first plot
    belief_plt = scatter!(
        left_panel,
        zeros(prey.observer.nBparticles),
        zeros(prey.observer.nBparticles),
        color=colour_posterior,
        strokewidth=0,
        markersize=size_belief,
    )
end  # plot internal particles

if PLOT_ARRAYS

    Likely_plt = plot!(
        middle_panel,
        OffsetArrays.no_offset_view(prey.observer.likelihood),
    )

    Posty_plt =
        plot!(right_panel, OffsetArrays.no_offset_view(prey.observer.posterior))

    predator_right_plt = poly!(
        right_panel,
        lift(
            s -> decompose(
                Point2f0,
                Circle(
                    Point2f0(
                        mat_radius + predator.x[],
                        mat_radius + predator.y[],
                    ),
                    predator.radius,
                ),
            ),
            t,
        ),
        color=RGBA(0.0, 0.0, 0.0, 0.0),
        strokecolor=predator.edgecolor,
        strokewidth=1.0,
    )
    predator_middle_plt = poly!(
        middle_panel,
        lift(
            s -> decompose(
                Point2f0,
                Circle(
                    Point2f0(
                        mat_radius + predator.x[],
                        mat_radius + predator.y[],
                    ),
                    predator.radius,
                ),
            ),
            t,
        ),
        color=RGBA(0.0, 0.0, 0.0, 0.0),
        strokecolor=predator.edgecolor,
        strokewidth=1.0,
    )

    prey_Lcopy_plt = poly!(
        middle_panel,
        decompose(
            Point2f0,
            Circle(Point2f0(mat_radius, mat_radius), prey.radius),
        ),
        color=prey.color,
        strokewidth=1,
        strokecolor=RGB(0.5, 0.75, 0.85),
    )
    prey_Pcopy_plt = poly!(
        right_panel,
        decompose(
            Point2f0,
            Circle(Point2f0(mat_radius, mat_radius), prey.radius),
        ),
        color=prey.color,
        strokewidth=1,
        strokecolor=RGB(0.5, 0.75, 0.85),
    )
end

# Prey
prey_plt = poly!(left_panel,
       decompose(Point2f0, Circle(Point2f0(0., 0.), prey.radius)),
       color=prey.color, strokewidth=1, strokecolor=RGB(.5, .75, .85))
preyGut_plt = poly!(left_panel,
      decompose(Point2f0, Circle(Point2f0(0., 0.), prey.gutradius)),
      color=prey.gutcolor, strokewidth=0.0)


receptor_plt = scatter!(left_panel, prey.receptor.x, prey.receptor.y ,
            markersize=prey.receptor.size,
            color=[prey.receptor.openColor for i in 1:prey.receptor.N],
            strokecolor=:black, strokewidth=0.25)

crystal_plt = scatter!(left_panel, prey.photoreceptor.x, prey.photoreceptor.y,
            markersize=prey.photoreceptor.size, marker=:diamond,
            color=[prey.photoreceptor.lightColor for i in 1:prey.photoreceptor.N],
            strokecolor=:black, strokewidth=0.25)


# reset axis limits (have been auto-adjusted by MakieLayout)
xlims!(left_panel, -mat_radius, mat_radius)
ylims!(left_panel, -mat_radius, mat_radius)
if PLOT_ARRAYS
    xlims!(middle_panel, 0, WorldSize)
    ylims!(middle_panel, 0, WorldSize)
    xlims!(right_panel, 0, WorldSize)
    ylims!(right_panel, 0, WorldSize)
end

ELECTRORECEPTION = true
PHOTORECEPTION = false

NTRIALS = [0]
tick()


FileName = "PlacozoanStalker" * string(Int(ELECTRORECEPTION)) * string(Int(PHOTORECEPTION)) * "_" *
    string(Nreceptors) * ".csv"

    CSV.write(FileName,
DataFrame(Range = Float64[], predatorx=Float64[], predatory = Float64[], xMAP = Int64[], yMAP = Int64[], 
               PosteriorEntropy=Float64[], KLD=Float64[], Nreceptors = Int[], 
               n_likelihood_particles=Int64[], n_posterior_particles = Int64[],  priorDensity = Float64[]))

for rep = 1:nReps
    for n_likelihood_particles = [1024 2048 4096 8192]
        for n_posterior_particles = n_likelihood_particles .÷ [1 2 4 8]
            for priorDensity = [0.0005 0.001 0.002 0.004 0.008]

                    # FileName = "PlacozoanStalker" * string(Int(ELECTRORECEPTION)) * string(Int(PHOTORECEPTION)) * "_" *
                    # string(Nreceptors) * "_" * string(n_likelihood_particles) * "_" * string(n_posterior_particles) *
                    # "_" * string(priorDensity) * "_" * string(rep)

        # record(scene, FileName * ".mp4", framerate=18, 1:nFrames) do i
                for i in 1:nFrames

    # predator random walk to within Δ of prey
                    stalk(predator, prey, approach_Δ)

    # electroreception
                    if ELECTRORECEPTION
                        electroreception(prey, predator)
        # set color of each receptor, indicating open or closed state
                        receptorColor = [prey.receptor.closedColor for j = 1:prey.receptor.N]
                        receptorColor[findall(x -> x == 1, prey.receptor.state)] .=
            prey.receptor.openColor
                        receptor_plt.color[] = receptorColor
                    end


    # photoreception
                    if PHOTORECEPTION
                        photoreception(prey, predator)
        # set color of each receptor, indicating open or closed state
                        crystalColor = [prey.photoreceptor.lightColor  for j in 1:prey.photoreceptor.N]
                        crystalColor[findall(x -> x == 1, prey.photoreceptor.state)] .=
            prey.photoreceptor.darkColor
                        crystal_plt.color[] = crystalColor
                    end


    # inference
    # ELECTRORECEPTION & PHOTORECEPTION are bools 
    #   specifying whether the respective receptor states are included in the inference
                    likelihood(prey, ELECTRORECEPTION, PHOTORECEPTION)  

                    sample_likelihood(prey)     # random sample from likelihood
                    bayesParticleUpdate(prey)   # Bayesian particle filter

                    (observation, belief) = reflect(prey) # reflect samples into margin

                    if PLOT_EXT_PARTICLES
        # update likelihood particle plot
                        Lparticle_plt[1] = prey.observer.Lparticle[:, 1]
                        Lparticle_plt[2] = prey.observer.Lparticle[:, 2]

        # update posterior particle plot
                        Bparticle_plt[1] = prey.observer.Bparticle[:, 1]
                        Bparticle_plt[2] = prey.observer.Bparticle[:, 2]
                    end

                    if PLOT_INT_PARTICLES

        # update observation particle plot
                        observation_plt[1] = observation[:, 1]
                        observation_plt[2] = observation[:, 2]

        # update observation particle plot
                        belief_plt[1] = belief[:, 1]
                        belief_plt[2] = belief[:, 2]

                    end

                    if PLOT_ARRAYS
                        bayesArrayUpdate(prey)
                        Likely_plt[1] =
            mask .* OffsetArrays.no_offset_view(prey.observer.likelihood)
                        Posty_plt[1] = mask .* OffsetArrays.no_offset_view(prey.observer.posterior)
                    end

    # record posterior entropy (& display during simulation)
                    prey.observer.PosteriorEntropy[i] = entropyBits(prey.observer)
                    prey.observer.KLD[i] = KLDBits(prey.observer)
                    recordRange(prey.observer, predator, i)
                    recordKLDBits(prey.observer, i)
    # println("t=", i, ", ", prey.observer.PosteriorEntropy[1] -
    #         prey.observer.PosteriorEntropy[i], ", ", prey.observer.KLD[i])


    # clock display
                    clock_plt.text =
        "                                                  t = " *
        string(floor(t[])) *
        "s"

    # Node update causes redraw
                    t[] = dt * (i + 1)


        # MAP predator location            
                    iMAP = findmax(prey.observer.posterior)[2]

# save data
                    CSV.write(FileName,
    DataFrame(hcat(prey.observer.range[i], predator.x[], predator.y[], iMAP[1], iMAP[2],
                   prey.observer.PosteriorEntropy[i], prey.observer.KLD[i], Nreceptors, 
                   n_likelihood_particles, n_posterior_particles,  priorDensity)),
                    header=false, append=true)

               print(".")


                end # frame
                println()

                NTRIALS[] = NTRIALS[] + 1
                laptimer()
                println(NTRIALS[], ", ", n_likelihood_particles, ", ", n_posterior_particles, ", ", priorDensity, ", ",  rep)
                initialize_particles(prey) # draw initial sample from prior
                initialize_prior(prey)     # initialize numerical Bayesian prior
                θ = π * rand()[] # Random initial heading (from above)
                predator.x[] = (mat_radius + 0.5 * predator_radius) * cos(θ)
                predator.y[] = (mat_radius + 0.5 * predator_radius) * sin(θ)
                t[] = 0

            end # n_prior_particles
        end # n_likelihood_particles
    end # priorDensity
end # rep


# end # main