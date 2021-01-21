# PlacozoanStalker
# function main(nFrames=100)

using TickTock
# import BayesianPlacozoan
# function main()

# to create a video file, find the line starting "record(scene ..."
# and comment it out (about line 400)
    
# choose what gets plotted (for demo/explanatory videos)
# the default is to plot 3 panels showing all particles + likelihood and posterior
# in which case the following switches are all ON (true)
PLOT_EXT_PARTICLES = true
PLOT_INT_PARTICLES = true
PLOT_ARRAYS = true

# DO_PLOTS switches plotting ON/OFF, for running multiple simulations
# to collect data without plotting. DO_PLOTS must be true for the 
# settings above to take effect
DO_PLOTS = true
if DO_PLOTS == false
    PLOT_EXT_PARTICLES = false
    PLOT_INT_PARTICLES = false
    PLOT_ARRAYS = false
end


# simulation parameters
nReps = 16
nFrames = 180       # number of animation frames
burn_time = 4         # compute initial prior by burning in with predator at "infinity"
mat_radius = 400
approach_Δ = 32.0         # predator closest approach distance
dt = 1.00

# # prey 
# priorDensity = 0.002
# posteriorSD = 100.0
# n_likelihood_particles = 5000
# n_posterior_particles = 2500


#  prey parameters
prey_radius = 120
prey_margin = 40
Nreceptors = 32
Ncrystals = 8
prey_fieldrange = 0   # no field


# predator parameters
predator_radius = 150
predator_margin = 0
predator_speed = 1.0
predator_fieldrange = mat_radius



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




# initialize prey observer

# likelihood(prey)           # initialize likelihood given initial receptor states
# sample_likelihood(prey)    # sample from normalized likelihood
# initialize_particles(prey) # draw initial sample from prior
# initialize_prior(prey)     # initialize numerical Bayesian prior

# time observable
# used to force scene update (nothing depends explicitly on time)
t = Node(0.0)

ELECTRORECEPTION = true
PHOTORECEPTION = false

NTRIALS = [0]



FileName = "PlacozoanStalker" * string(Int(ELECTRORECEPTION)) * string(Int(PHOTORECEPTION)) * "_" *
    string(Nreceptors) 

CSV.write(FileName * ".csv",
    DataFrame(Range=Float64[], predatorx=Float64[], predatory=Float64[], xMAP=Int64[], yMAP=Int64[], 
                   PosteriorEntropy=Float64[], KLD=Float64[], 
                   Dmed = Float64[], Dquart = Float64[], D5pc = Float64[], D1pc = Float64[],
                   Dmin = Float64[], 
                   Θ1pc = Float64[], Θ5pc = Float64[], Θquart = Float64[], Θmed = Float64[],
                   Θquarta = Float64[], Θ5pca = Float64[], Θ1pca = Float64[],     
                   Θmin = Float64[], Θmax = Float64[],           
                   Nreceptors=Int[], n_likelihood_particles=Int64[], n_posterior_particles=Int64[],  
                   priorDensity=Float64[]))

# use dummy predator and prey to precompute fields, receptive fields and prior
# these are common to all Placozoans with the same number of receptors
# so only need to be computed once
dummy_prey = Placozoan(prey_radius, prey_margin, prey_fieldrange,
    Nreceptors, sizeof_receptor, mat_radius,
    Ncrystals, sizeof_crystal, mat_radius,
    1, 1, 1, 1)



dummy_predator = Placozoan(predator_radius, predator_margin, predator_fieldrange,
RGBA(.25, 0.1, 0.1, 1.0),
RGBA(.45, 0.1, 0.1, 0.25),
RGB(.95, 0.1, 0.1) )
dummy_predator.x[] = 0
dummy_predator.y[] = mat_radius + 2.0*predator_radius
dummy_predator.speed[] = predator_speed

placozoanFieldstrength!(dummy_predator)
Ereceptor_RF(dummy_prey, dummy_predator)
Vreceptor_RF(dummy_prey)
initialize_prior(dummy_prey)


# burn in
for i in 1:burn_time
    if ELECTRORECEPTION
        electroreception(dummy_prey, dummy_predator)
    end
    if PHOTORECEPTION
        photoreception(dummy_prey, dummy_predator)
    end
    likelihood(dummy_prey, ELECTRORECEPTION, PHOTORECEPTION)  
    bayesArrayUpdate(dummy_prey)
    dummy_prey.observer.prior[:,:] = dummy_prey.observer.posterior[:,:]
end

radialSmooth(dummy_prey.observer.prior, prey_radius:mat_radius)
# # burn in
# for i in 1:dummy_prey.observer.burnIn
#     if ELECTRORECEPTION
#         electroreception(dummy_prey, dummy_predator)
#     end
#     if PHOTORECEPTION
#         photoreception(dummy_prey, dummy_predator)
#     end
#     likelihood(dummy_prey, ELECTRORECEPTION, PHOTORECEPTION)  
#     bayesArrayUpdate(dummy_prey)
#     dummy_prey.observer.prior[:,:] = (1.0 - 1/i)*dummy_prey.observer.prior[:,:] + dummy_prey.observer.posterior[:,:]/i
# end




tick()

for rep = 1 # 1:nReps
    for n_likelihood_particles = [2048] # [512 1024  2048 4096 8192 ]
        for n_posterior_particles = [1024] # Int.(n_likelihood_particles .÷ [.5 1 2 4])
            for posteriorDeaths = [32] # [.001 .01 .1]

                # construct placozoans
                # HINT: These are local variables but if they are declared global 
                # then they appear in the REPL workspace if the program is interrupted (Ctrl-C)
                global prey = Placozoan(prey_radius, prey_margin, prey_fieldrange,
                    Nreceptors, sizeof_receptor, mat_radius,
                    Ncrystals, sizeof_crystal, mat_radius,
                    n_likelihood_particles, Int(n_posterior_particles),
                    posteriorDeaths, nFrames)

                prey.receptor.pOpen[:] = dummy_prey.receptor.pOpen[:] 
                prey.photoreceptor.pOpenV[:] = dummy_prey.photoreceptor.pOpenV[:] 
                # println("____")
                # println(maximum(prey.receptor.pOpen[1]), " ", maximum(dummy_prey.receptor.pOpen[1]))
                # println("____")
                # if NTRIALS[]==0
               # initialize_prior(prey)     # initialize numerical Bayesian prior
               prey.observer.prior[:,:] = dummy_prey.observer.prior[:,:]
               prey.observer.posterior[:,:] = prey.observer.prior[:,:]
               initialize_particles(prey) # draw initial sample from prior
               # end
                
                # initializeObserver(prey, n_likelihood_particles, n_posterior_particles, priorDensity)

                # predator constructed without observer
                global predator = Placozoan(predator_radius, predator_margin, predator_fieldrange,
                     RGBA(.25, 0.1, 0.1, 1.0),
                     RGBA(.45, 0.1, 0.1, 0.25),
                     RGB(.95, 0.1, 0.1) )
                predator.speed[] = predator_speed
                θ = π * rand()[] # Random initial heading (from above)
                predator.x[] = (mat_radius + 0.5 * predator_radius) * cos(θ)
                predator.y[] = (mat_radius + 0.5 * predator_radius) * sin(θ)
                predator.field[:] = dummy_predator.field[:]
                predator.potential[:] = dummy_predator.potential[:]




                if DO_PLOTS
                
                    if !PLOT_ARRAYS    # not plotting likelihoods or posterior
                        scene, layout = layoutscene(resolution=(WorldSize, WorldSize))
                        left_panel =    layout[1, 1] = LAxis(scene, title="Placozoan", 
                                                         backgroundcolor=colour_background )
                # clock_panel = LAxis(scene,
                #      bbox = BBox(100 , 100, 100, 100) )
                        hidespines!(left_panel)
                        hidedecorations!(left_panel)
                    else
                        scene, layout = layoutscene(resolution = ( Int(round(3 * .75 * WorldSize)), Int(round(.75 * WorldSize + 40)) ) )
                        shim = layout[1, 1] = LAxis(scene)
                        left_panel = layout[1, 2] = 
                            LAxis( scene,  title="Placozoan: ( " * string(n_likelihood_particles) * ", " *
                                string(n_posterior_particles) * ", " * string(posteriorDeaths) * " )",
                                backgroundcolor=colour_background )
                        middle_panel = layout[1, 3] = LAxis(scene, title="Likelihood")
                        right_panel = layout[1, 4] = LAxis(scene, title="Bayesian Observer")
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
                    display(scene)
                end # DO_PLOTS
                
                if PLOT_ARRAYS
                    mat_middle_plt = poly!(middle_panel,
                        decompose(Point2f0, Circle(Point2f0(mat_radius, mat_radius), mat_radius)),
                        color=RGBA(0.0, 0.0, 0.0, 0.0),
                        strokewidth=0.5,
                        strokecolor=RGB(0.75, 0.75, 0.75),
                    )

                    mat_right_plt = poly!(right_panel,
                        decompose(Point2f0, Circle(Point2f0(mat_radius, mat_radius), mat_radius)),
                        color=RGBA(0.0, 0.0, 0.0, 0.0),
                        strokewidth=0.5,
                        strokecolor=RGB(0.75, 0.75, 0.75),
                    )
                end  # PLOT_ARRAYS
                
                if DO_PLOTS
                    # display nominal time on background
                    clock_plt = LText(scene, "                         t = 0.0s",
                        color=:white, textsize=18, halign=:left)

                # predator drawn using lift(..., node)
                # (predatorLocation does not depend explicitly on t, but this causes
                #  the plot to be updated when the node t changes)
                predator_plt = poly!(left_panel,
                    lift(s -> decompose(Point2f0, Circle(Point2f0(predator.x[], predator.y[]),
                        predator.radius)), t),
                    color=predator.color, strokecolor=predator.edgecolor, strokewidth=.5)

                end # DO_PLOTS


                if PLOT_EXT_PARTICLES

                    # plot likelihood particles (samples from likelihood)
                    Lparticle_plt = scatter!(left_panel,
                        prey.observer.Lparticle[1:prey.observer.nLparticles[], 1],
                        prey.observer.Lparticle[1:prey.observer.nLparticles[], 2],
                        color=:yellow, markersize=size_likelihood, strokewidth=0.1)

                    # plot posterior particles
                    Bparticle_plt = scatter!(left_panel,
                        prey.observer.Bparticle[1:prey.observer.nBparticles[], 1],
                        prey.observer.Bparticle[1:prey.observer.nBparticles[], 2],
                        color=colour_posterior, markersize=size_posterior, strokewidth=0.1)

                end # plot external particles


                if PLOT_INT_PARTICLES

                # plot projection of likelihood particles into prey margin
                # nb this is a dummy plot
                # the correct particle locations are inserted before first plot
                observation_plt = scatter!(left_panel,
                    zeros(prey.observer.nLparticles[]),
                    zeros(prey.observer.nLparticles[]),
                    color=:yellow, strokewidth=0, markersize=size_observation )

                # plot projection of posterior particles into prey margin
                # nb this is a dummy plot
                # the correct particle locations are inserted before first plot
                belief_plt = scatter!(left_panel,
                    zeros(prey.observer.nBparticles[]),
                    zeros(prey.observer.nBparticles[]),
                    color=colour_posterior, strokewidth=0, markersize=size_belief)

                end  # plot internal particles

                if PLOT_ARRAYS

                    Likely_plt = plot!( middle_panel,
                        OffsetArrays.no_offset_view(prey.observer.likelihood), colorrange = (0.0, 1.0), colormap = :haline)

                    Posty_plt =  surface!(right_panel, 1:WorldSize, 1:WorldSize,
                        OffsetArrays.no_offset_view(prey.observer.posterior),  colormap = :plasma)

                    PostContour_plt = contour!(right_panel, 1:WorldSize, 1:WorldSize,
                        lift(u->u, Posty_plt[3]), levels = [1.0e-6, 1.0e-5, 1.0e-4], color = RGB(.35,.35,.35))

                    predator_right_plt = poly!(right_panel,
                        lift(s -> decompose(Point2f0, Circle(Point2f0(
                            mat_radius + predator.x[], mat_radius + predator.y[]),
                            predator.radius)), t ),
                            color=RGBA(0.0, 0.0, 0.0, 0.0),
                            strokecolor=predator.edgecolor,
                            strokewidth=1.0,
                        )

                    predator_middle_plt = poly!( middle_panel,
                        lift(s -> decompose(Point2f0,
                            Circle(Point2f0(mat_radius + predator.x[],mat_radius + predator.y[]),
                            predator.radius) ), t),
                            color=RGBA(0.0, 0.0, 0.0, 0.0),
                            strokecolor=predator.edgecolor,
                            strokewidth=1.0    
                        )

                    prey_Lcopy_plt = poly!(middle_panel,
                        decompose(Point2f0, Circle(Point2f0(mat_radius, mat_radius), prey.radius)),
                        color=:black, strokewidth=1, strokecolor=RGB(0.0, 0.0, 0.0))

                    prey_Pcopy_plt = poly!(right_panel,
                        decompose(Point2f0, Circle(Point2f0(mat_radius, mat_radius), prey.radius)),
                        color=:black, strokewidth=1, strokecolor=RGB(0.0, 0.0, 0.0))

                end  # plot arrays

                if DO_PLOTS

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

                    L_receptor_plt = scatter!(middle_panel, 
                        mat_radius .+1 .+prey.receptor.x, mat_radius .+1 .+prey.receptor.y ,
                        markersize=prey.receptor.size,
                        color=[prey.receptor.openColor for i in 1:prey.receptor.N],
                        strokecolor=:black, strokewidth=0.25)

                    R_receptor_plt = scatter!(right_panel, 
                        mat_radius .+1 .+prey.receptor.x, mat_radius .+1 .+prey.receptor.y ,
                        markersize=prey.receptor.size,
                        color=[prey.receptor.openColor for i in 1:prey.receptor.N],
                        strokecolor=:black, strokewidth=0.25)
                        
                    if PHOTORECEPTION
                    crystal_plt = scatter!(left_panel, prey.photoreceptor.x, prey.photoreceptor.y,
                        markersize=prey.photoreceptor.size, marker=:diamond,
                            color=[prey.photoreceptor.lightColor for i in 1:prey.photoreceptor.N],
                            strokecolor=:black, strokewidth=0.25)
                    end


                    # reset axis limits (have been auto-adjusted by MakieLayout)
                    xlims!(left_panel, -mat_radius, mat_radius)
                    ylims!(left_panel, -mat_radius, mat_radius)
                    if PLOT_ARRAYS
                        xlims!(middle_panel, 0, WorldSize)
                        ylims!(middle_panel, 0, WorldSize)
                        xlims!(right_panel, 0, WorldSize)
                        ylims!(right_panel, 0, WorldSize)
                    end

                end #DO_PLOTS

                videoName = FileName * "_" * string(rep) * "_" * string(n_likelihood_particles) * "_" *
                           string(n_posterior_particles) * "_" * string(posteriorDeaths) * ".mp4"

                # VIDEO RECORDING
                # comment out ONE of the following 2 lines to (not) generate video file
                record(scene, videoName , framerate=9, 1:nFrames) do i     # generate video file
                #for i in 1:nFrames                                      # just compute

                   #println(i)

                    # predator random walk to within Δ of prey
                    stalk(predator, prey, approach_Δ)

                    # electroreception
                    if ELECTRORECEPTION
                        electroreception(prey, predator)

                        if DO_PLOTS
                            # set color of each receptor, indicating open or closed state
                            receptorColor = [prey.receptor.closedColor for j = 1:prey.receptor.N]
                            receptorColor[findall(x -> x == 1, prey.receptor.state)] .= prey.receptor.openColor
                            receptor_plt.color[] = L_receptor_plt.color[] = R_receptor_plt.color[] = receptorColor
                        end # DO_PLOTS
                    end  # electroreception


                    # photoreception
                    if PHOTORECEPTION
                        photoreception(prey, predator)
                        if DO_PLOTS
                            # set color of each receptor, indicating open or closed state
                            crystalColor = [prey.photoreceptor.lightColor  for j in 1:prey.photoreceptor.N]
                            crystalColor[findall(x -> x == 1, prey.photoreceptor.state)] .= prey.photoreceptor.darkColor
                            crystal_plt.color[] = crystalColor
                        end # DO_PLOTS
                    end # photoreception


                    # inference
                    # ELECTRORECEPTION & PHOTORECEPTION are bools 
                    #   specifying whether the respective receptor states are included in the inference

                    likelihood(prey, ELECTRORECEPTION, PHOTORECEPTION)  

                    sample_likelihood(prey)     # random sample from likelihood
    
                    bayesParticleUpdate(prey)   # Bayesian particle filter

                    bayesArrayUpdate(prey)  # numerical sequential Bayes (benchmark)

                    (observation, belief) = reflect(prey) # reflect samples into margin
    

                    if PLOT_EXT_PARTICLES
                        # update likelihood particle plot
                        Lparticle_plt[1] = prey.observer.Lparticle[1:prey.observer.nLparticles[], 1]
                        Lparticle_plt[2] = prey.observer.Lparticle[1:prey.observer.nLparticles[], 2]

                        # update posterior particle plot
                        Bparticle_plt[1] = prey.observer.Bparticle[1:prey.observer.nBparticles[], 1]
                        Bparticle_plt[2] = prey.observer.Bparticle[1:prey.observer.nBparticles[], 2]
                    end # PLOT_EXT_PARTICLES

                    if PLOT_INT_PARTICLES

                        # update observation particle plot
                        observation_plt[1] = observation[1:prey.observer.nLparticles[], 1]
                        observation_plt[2] = observation[1:prey.observer.nLparticles[], 2]

                        # update observation particle plot
                        belief_plt[1] = belief[1:prey.observer.nBparticles[], 1]
                        belief_plt[2] = belief[1:prey.observer.nBparticles[], 2]

                    end # PLOT_INT_PARTICLES

                    if PLOT_ARRAYS
                        

                        Likely_plt[1] = mask .* OffsetArrays.no_offset_view(prey.observer.likelihood)
                        Posty_plt[3] = mask .* OffsetArrays.no_offset_view(prey.observer.posterior)
                    end # PLOT_ARRAYS

                    # record posterior entropy (& display during simulation)
                    prey.observer.PosteriorEntropy[i] = entropyBits(prey.observer)
                    prey.observer.KLD[i] = KLDBits(prey.observer)
                    recordRange(prey.observer, predator, i)
                    recordKLDBits(prey.observer, i)

                    # clock display
                    if DO_PLOTS
                        clock_plt.text =  "                         t = " *   
                                        string(Int(floor(t[]+1))) *  "s"

                        # Node update causes redraw
                        t[] = dt * (i + 1)

                    end # DO_PLOTS

                    # MAP predator location            
                    iMAP = findmax(prey.observer.posterior)[2]
                   # println(iMAP[1], ", ", iMAP[2])

                    (QD, Dmin, Qθ, θmin, θmax) = particleStats(prey, atan(predator.y[],predator.x[]) )
                    # println(QD, ", ", Dmin, ", ", Qθ, ", ", θmin, ", ", θmax)
                    # sleep(2)

                    # save data
                    CSV.write(FileName * ".csv",
                        DataFrame(hcat(prey.observer.range[i], predator.x[], predator.y[], iMAP[1], iMAP[2],
                            prey.observer.PosteriorEntropy[i], prey.observer.KLD[i], 
                            QD..., Dmin, Qθ..., θmin,  θmax, 
                            Nreceptors, 
                            n_likelihood_particles, n_posterior_particles,  posteriorDeaths)),
                            header=false, append=true)

                    print(".")


                end # frame

                println()

                NTRIALS[] = NTRIALS[] + 1
                laptimer()
                println(NTRIALS[], ", ", n_likelihood_particles, ", ", n_posterior_particles, ", ", posteriorDeaths, ", ",  rep)
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









   
