# PlacozoanStalker
# function main(nFrames=100)

using TickTock
# import BayesianPlacozoan
# function main()

# to [not] create a video file, uncomment the line starting "record(scene ..."
# (about line 400) and comment out the for i = ... (next line)
    
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
nReps = 64
nFrames = 480       # number of animation frames
burn_time = 30      # burn in posterior initially for 30 sec with predator outside observable world
mat_radius = 400
approach_Δ = 25.0         # predator closest approach distance
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
predator_speed = 0.25
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


# open file to save results as dataframe, 1 row per simulation step (=per video frame when videoed)
FileName = "PlacozoanStalker" * string(Int(ELECTRORECEPTION)) * string(Int(PHOTORECEPTION)) * "_" *
    string(Nreceptors) 

CSV.write(FileName * ".csv",
    DataFrame(rep=Int64[],
    
        # distance to closest edge of predator, and predator location
        Range=Float64[], predatorx=Float64[], predatory=Float64[], 
        
        # Bayesian MAP estimate of predator location
        xMAP=Int64[], yMAP=Int64[], 

        # Entropy of Bayesian posterior & Kullback-Leibler divergence (relative entropy)
        #   of particle distribution (KLD), a uniform random sample (KLD0) and a sample 
        # from the posterior density (KLDI, ie simulating an ideal Bayesian particle filter),
        # with respect to the posterior (ie measure information loss in the particle estimate,
        #   a sample from the posterior and a random sample, relative to the posterior)
        PosteriorEntropy=Float64[], KLD=Float64[], KLD0=Float64[], KLDI=Float64[],

        ## posterior density summary stats ##

        # posterior probability that predator is closer than 25, 50 and 100um
        PR25 = Float64[], PR50 = Float64[], PR100 = Float64[], 

        # quantiles of posterior density of distance to predator
        # 1%, 5%, 25% and 50% (median) (proximal side only, ie we care about how close the predator
        # might be, not how far away it might be)
        QP01 = Float64[], QP05 = Float64[], QP25 = Float64[], QP50 = Float64[],

        # quantiles of angular distribution
        Qψ01 = Float64[], Qψ05 = Float64[], Qψ25 = Float64[], Qψ50 = Float64[], 
        Qψ75 = Float64[], Qψ95 = Float64[], Qψ99 = Float64[],

        ## Particle summary stats  ##

        # proportion of particles within 25, 50 and 100um
        NR25 = Float64[], NR50 = Float64[], NR100 = Float64[], 

        # quantiles of particle proximity, 1%, 5%, 25% and 50%
        # e.g. QN05 is range including closest 5% of particles
        QN01 = Float64[], QN05 = Float64[], QN25 = Float64[], QN50 = Float64[], 

        # quantiles of particle direction error (from true heading to predator)
        # giving credibility intervals for direction
        # e.g. QΘ05 is left/anticlockwise limit of 5% credibility interval for direction to predator,
        #      and Q095 is right/clockwise limit
        QΘ01 = Float64[], QΘ05 = Float64[], QΘ25 = Float64[], QΘ50 = Float64[], 
        QΘ75 = Float64[], QΘ95 = Float64[], QΘ99 = Float64[],

        # M-cell posterior belief of predator in patch
        MP = Float64[],
        
        # trial parameters
        Nreceptors=Int[], n_likelihood_particles=Int64[], n_posterior_particles=Int64[],  
        priorDensity=Float64[]))


tick()

for rep = 1:nReps
    for n_likelihood_particles = [1024] # [512 1024  2048 4096 8192 ]
        for n_posterior_particles = [1024] # Int.(n_likelihood_particles .÷ [.5 1 2 4])
            for posteriorDeaths = [4] # [.001 .01 .1]

                # construct placozoans
                # HINT: These are local variables but if they are declared global 
                # then they appear in the REPL workspace if the program is interrupted (Ctrl-C)
                global
                prey = Placozoan(prey_radius, prey_margin, prey_fieldrange,
                    Nreceptors, sizeof_receptor, mat_radius,
                    Ncrystals, sizeof_crystal, mat_radius,
                    n_likelihood_particles, Int(n_posterior_particles),
                    posteriorDeaths, nFrames)

                # predator constructed without observer
                global
                predator = Placozoan(predator_radius, predator_margin, predator_fieldrange,
                     RGBA(.25, 0.1, 0.1, 1.0),
                     RGBA(.45, 0.1, 0.1, 0.25),
                     RGB(.95, 0.1, 0.1) )
                predator.speed[] = predator_speed
                θ = π * rand()[] # Random initial heading (from above)
                predator.x[] = (mat_radius + predator_radius) * cos(θ)
                predator.y[] = (mat_radius + predator_radius) * sin(θ)
                # predator.field[:] = dummy_predator.field[:]
                # predator.potential[:] = dummy_predator.potential[:]


            
            placozoanFieldstrength!(predator)
            Ereceptor_RF(prey, predator)
            Vreceptor_RF(prey)
            initialize_prior(prey)
            prey.observer.posterior[:,:] = prey.observer.prior[:,:]

            # burn in posterior
            # from uniform to posterior given no predator in the observable world for burn-in time
            for i in 1:burn_time
                if ELECTRORECEPTION
                    electroreception(prey, predator)
                end
                if PHOTORECEPTION
                    photoreception(prey, predator)
                end
                likelihood(prey, ELECTRORECEPTION, PHOTORECEPTION)  
                bayesArrayUpdate(prey)
                prey.observer.prior[:,:] = prey.observer.posterior[:,:]
                radialSmooth(prey.observer.prior, prey_radius:mat_radius)           
            end
      
            predator.x[] = (mat_radius + 0.0* predator_radius) * cos(θ)
            predator.y[] = (mat_radius + 0.0 * predator_radius) * sin(θ)                



            # initialize particle filter
            initialize_particles(prey) # draw initial sample from prior

                if DO_PLOTS
                
                    if !PLOT_ARRAYS    # not plotting likelihoods or posterior
                        scene, layout = layoutscene(resolution=(WorldSize, WorldSize))
                        left_panel =    layout[1, 1] = LAxis(scene, title="Placozoan", 
                                                         backgroundcolor=colour_background )

                        hidespines!(left_panel)
                        hidedecorations!(left_panel)
                    else
                        scene = Figure(resolution = ( Int(round(3 * .75 * WorldSize)), Int(round(.75 * WorldSize + 40)) ), backgroundcolor = :black)
                        left_panel = scene[1,1] = Axis(scene, title="Placozoan Model (Bayesian Particle Filter)",
                            titlecolor = title_color, backgroundcolor=:black )
                        middle_panel = scene[1, 2] = Axis(scene, title="Likelihood", titlecolor = title_color, backgroundcolor=colour_background)
                        right_panel = scene[1, 3] = Axis(scene, title="Bayesian Observer", titlecolor = title_color, backgroundcolor=colour_background)


                        
                        # colsize!(scene, 1, Relative(0.04))
                        # colsize!(scene, 2, Relative(0.32))
                        # colsize!(scene, 3, Relative(0.32))
                        # colsize!(scene, 4, Relative(0.32))

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
                    Δ_plt = Label(scene, "               Δ = 0",
                        color=RGB(.7,.7, .7), textsize=18, halign=:left)


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
                        color=colour_likelihood, markersize=size_likelihood, strokewidth=0.1)

                    # plot posterior particles
                    Pparticle_plt = scatter!(left_panel,
                        prey.observer.Pparticle[1:prey.observer.nPparticles[], 1],
                        prey.observer.Pparticle[1:prey.observer.nPparticles[], 2],
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
                    zeros(prey.observer.nPparticles[]),
                    zeros(prey.observer.nPparticles[]),
                    color=colour_posterior, strokewidth=0, markersize=size_belief)

                end  # plot internal particles

                if PLOT_ARRAYS

                    Likely_plt = plot!( middle_panel,
                        OffsetArrays.no_offset_view(prey.observer.likelihood), colorrange = (0.0, 1.25), colormap = :copper)  # :turku

                    Posty_plt =  surface!(right_panel, 1:WorldSize, 1:WorldSize,
                        OffsetArrays.no_offset_view(prey.observer.posterior),  colormap = :magma)

                    # PostContour_plt = contour!(right_panel, 1:WorldSize, 1:WorldSize,
                    #     lift(u->u, Posty_plt[3]), levels = [1.0e-6, 1.0e-5, 1.0e-4], color = RGB(.35,.35,.35))

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
                        color=RGBA(0.0, 0.0, 0.0, 0.0), strokewidth=1.5, strokecolor=prey.gutcolor)

                    preygut_Lcopy_plt = poly!(middle_panel,
                        decompose(Point2f0, Circle(Point2f0(mat_radius, mat_radius), prey.radius - prey.marginwidth + 1)),
                        color=prey.gutcolor, strokewidth=0.5, strokecolor=prey.gutcolor*.75)

                    prey_Pcopy_plt = poly!(right_panel,
                        decompose(Point2f0, Circle(Point2f0(mat_radius, mat_radius), prey.radius)),
                        color=RGBA(0.0, 0.0, 0.0, 0.0), strokewidth=1.5, strokecolor=prey.gutcolor)

                    preygut_Pcopy_plt = poly!(right_panel,
                        decompose(Point2f0, Circle(Point2f0(mat_radius, mat_radius), prey.radius - prey.marginwidth + 1)),
                        color=prey.gutcolor, strokewidth=0.5, strokecolor=prey.gutcolor*.75)


                end  # plot arrays

                if DO_PLOTS

                    # Prey
                    prey_plt = poly!(left_panel,
                        decompose(Point2f0, Circle(Point2f0(0., 0.), prey.radius)),
                        color=prey.color, strokewidth=0.25, strokecolor=prey.gutcolor)

                    preyGut_plt = poly!(left_panel,
                        decompose(Point2f0, Circle(Point2f0(0., 0.), prey.gutradius)),
                        color=prey.gutcolor, strokewidth=0.5, strokecolor=prey.gutcolor*.75)


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

                videoName = FileName * "_" * string(n_likelihood_particles) * "_" *
                           string(n_posterior_particles) * "_" * string(posteriorDeaths) * "_" * string(rep) * ".mp4"

                # VIDEO RECORDING
                # comment out ONE of the following 2 lines to (not) generate video file
                #record(scene, videoName , framerate=16, 1:nFrames) do i     # generate video file
                for i in 1:nFrames                                      # just compute

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

                    #(observation, belief) = reflect!(prey) # reflect samples into margin
    

                    if PLOT_EXT_PARTICLES
                        # update likelihood particle plot
                        Lparticle_plt[1] = prey.observer.Lparticle[1:prey.observer.nLparticles[], 1]
                        Lparticle_plt[2] = prey.observer.Lparticle[1:prey.observer.nLparticles[], 2]

                        # update posterior particle plot
                        Pparticle_plt[1] = prey.observer.Pparticle[1:prey.observer.nPparticles[], 1]
                        Pparticle_plt[2] = prey.observer.Pparticle[1:prey.observer.nPparticles[], 2]
                    end # PLOT_EXT_PARTICLES

                    if PLOT_INT_PARTICLES

                        # update observation particle plot
                        observation_plt[1] = prey.observer.Sparticle[1:prey.observer.nLparticles[], 1]
                        observation_plt[2] = prey.observer.Sparticle[1:prey.observer.nLparticles[], 2]

                        # update observation particle plot
                        belief_plt[1] = prey.observer.Bparticle[1:prey.observer.nPparticles[], 1]
                        belief_plt[2] = prey.observer.Bparticle[1:prey.observer.nPparticles[], 2]

                    end # PLOT_INT_PARTICLES

                    if PLOT_ARRAYS
                        

                        Likely_plt[1] = mask .* OffsetArrays.no_offset_view(prey.observer.likelihood)
                        Posty_plt[3] = mask .* OffsetArrays.no_offset_view(prey.observer.posterior)
                    end # PLOT_ARRAYS

                    # record posterior entropy (& display during simulation)
                    prey.observer.PosteriorEntropy[i] = entropy(prey.observer.posterior)
                    #Sprey.observer.KLD[i] = KLDBits(prey.observer)
                    recordRange(prey.observer, predator, i)
                    KLD!(prey.observer, i)

                    # clock display
                    if DO_PLOTS
                        Δ_plt.text =  "               Δ = " *  
                          string(Int(round(sqrt(predator.x[]^2+predator.y[]^2)-predator.radius-prey.radius)))

                        # Node update causes redraw
                        t[] = dt * i

                    end # DO_PLOTS

                    # MAP predator location            
                    iMAP = findmax(prey.observer.posterior)[2]
                   # println(iMAP[1], ", ", iMAP[2])

                    # summary stats of particle distribution
                    # NR = proportion of particles (estimated probabilty) that predator is closer than 25,50 & 100um
                    # QN = quantiles of particle range, [0.005 0.025 0.25 0.5 0.75 0.975 0.995]
                    #      giving 1%, 5% and 50% credibility intervals + median estimate
                    # QΘ = quantiles of particle angle deviation from heading to predator (as above)
                    # MP = M-cell's posterior belief that there is a predator in its patch 
                    (NR, QN, Qθ, MP) = particleStats(prey, predator) 

                    (PR, QP, Qψ) = observerStats(prey, predator) 

                 
                    # println(QD, ", ", Dmin, ", ", Qθ, ", ", θmin, ", ", θmax)
                    # sleep(2)

                    # save data (see file open command for more detailed description of variables saved)
                    CSV.write(FileName * ".csv",
                        DataFrame(hcat(rep, 
                        
                        # range and location 
                        prey.observer.range[i], predator.x[], predator.y[], iMAP[1], iMAP[2],

                        # entropy/information in particle filter
                        prey.observer.PosteriorEntropy[i], prey.observer.KLD[i], prey.observer.KLD0[i], prey.observer.KLDI[i], 
                        
                        # summary stats of posterior probability 
                        PR..., QP, Qψ, 

                        # summary stats of particle distribution
                        NR..., QN..., Qθ..., MP, 

                        # trial parameters
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
                predator.x[] = (mat_radius + 0.75 * predator_radius) * cos(θ)
                predator.y[] = (mat_radius + 0.75 * predator_radius) * sin(θ)
                t[] = 0
            
            end # n_prior_particles

        
        end # n_likelihood_particles
    
    end # priorDensity

end # rep


# end # main









   
