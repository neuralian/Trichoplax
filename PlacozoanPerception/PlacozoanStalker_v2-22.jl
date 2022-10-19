# PlacozoanStalker v1.0
# This file is essentially a script that calls functions from BayesianPlacozoan.jl
# and handles the display and data logging over multiple trials.
#
# parameters (specified in this file): 
#    N_LIKELIHOOD_PARTICLES: likelihood (sensory) particle population sizes (1d array)
#    N_POSTERIOR_PARTICLES:  posterior (belief) particle population sizes (1d array)
#    POSTERIOR_DEATH_RATE:   proportion of posterior particles replaced by samples from 
#                            the initial prior ('evolved bias') distribution on each iteration.
#                            This helps to prevent the posterior spuriuosly condensing 
#                            (hallucination)  and maintains 'vigilance' on the horizon
#                            when there is limited sensory information. Should be small 
#                            numbers e.g. [.0001] for 1 in every 1000 (1d array).
#    N_REPS:                 Number of replicate trials for each of the ebove parameters.
#    SHOW_ANIMATION:         true/false.
#                            NB: To record (or not) animations as .mpg video files, 
#                            find  "VIDEO RECORDING ON/OFF"  (in a comment about line 460) 
#                            and follow instructions there.  SHOW_ANIMATION must be true
#                            to record animations (ie there must be animations to record).
#    LOG_DATA:               Log simulation parameters and data to .csv file.
#                            (Filenames are created from date/time & parameters).
#  NB: When collecting simulation data looping over multiple parameter values it can be useful to show 
#      animations as a way to check that the simulation is running OK. This slows down the simulation
#      a bit, but not so much in the grand scheme & IMO verification is more important than speed. However,
#      it is probably quite a bad idea to record all of these animations; they are large files 
#      and it may take you hours to look at them!  Better to create videos for particular parameter values
#      (based on data analysis indicating 'interesting' cases).
#
# There are various other parameters that can be tweaked, between here and the start of the 
# simulation loop.
#
#  (c) Michael G Paulin neuralian@gmail.com 2020-2022  


using TickTock
using Dates
include("BayesianPlacozoan_v2-22.jl")

# wrap the script in a function so it compiles to (more) specialized (faster) code
function placozoanStalker()


# SIMULATION parameters
n_likelihood_particles =  6400
n_posterior_particles = 12800
posteriorDeathRate = .0005

# number of replicate simulations for each combination of the above parameters
N_REPS = 32

# show animation while simulating true/false
# (must be true if )
SHOW_ANIMATION = true

# log parameters and simulation statistics per trial true/false 


LOG_DATA = false



# choose what gets plotted (for demo/explanatory videos)
# the default is to plot 3 panels showing all particles + likelihood and posterior
# in which case the following switches are all ON (true)
PLOT_EXT_PARTICLES = true
PLOT_INT_PARTICLES = true
PLOT_ARRAYS = true


if SHOW_ANIMATION == false
    PLOT_EXT_PARTICLES = false
    PLOT_INT_PARTICLES = false
    PLOT_ARRAYS = false
end


# simulation parameters
nFrames = 480    # number of animation frames
burn_time = 30      # burn in posterior initially for 30 sec with predator outside observable world
mat_radius = 400
min_Δ = 5.0         # predator closest approach distance
Δinit = 175.0
dt = 1.00

# # prey 
# priorDensity = 0.002
# posteriorSD = 100.0
# n_likelihood_particles = 5000
# n_posterior_particles = 2500


#  prey parameters-
prey_radius = 120
prey_margin = 40
Nreceptors = 32
Ncrystals = 32
prey_fieldrange = 0   # no field
prey_spinrate =   0.0005π  # prey turns slowly clockwise
particle_collisionRadius = 1.5 # when particles are this close to each other they are deemed to have collided
posterior_diffusion_coefficient = 4.0


# predator parameters
predator_radius = 150
predator_margin = 0
predator_speed =  .365
predator_brownian = 0.0 # 0.5*posterior_diffusion_coefficient
predator_fieldrange = mat_radius

# "Mauthner" cell

mcell_inset = 2.0*mcell_radius  # inset of M-cell centre from edge of animal



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


# time observable
# used to force scene update (nothing depends explicitly on time)
t = Observable{Int}(1)


ELECTRORECEPTION = true
PHOTORECEPTION = false

trialNumber = 0

# filename encodes number of receptors and number of frames (simulation timesteps) 
DataFileName = "PlacozoanStalkerV2.0 " * Dates.format(Dates.now(), "yyyy-mm-dd HH.MM")

# open file to save results as dataframe, 1 row per simulation step (=per video frame)
if LOG_DATA

    CSV.write(DataFileName * ".csv",
        DataFrame(rep=Int64[],

            # trial parameters
            n_likelihood_particles = Int64[], 
            n_posterior_particles  = Int64[], 
            n_vigilant_particles   = Float64[],
        
            # distance to closest edge of predator, and predator location
            Range     = Float64[], 
            predatorx = Float64[], 
            predatory = Float64[], 
            
            # Bayesian MAP estimate of predator location
            xMAP=Int64[], 
            yMAP=Int64[], 

            # Entropy of Bayesian posterior & Kullback-Leibler divergence (relative entropy)
            #   of particle distribution (KLD), a uniform random sample (KLD0) and a sample 
            # from the posterior density (KLDI, ie simulating an ideal Bayesian particle filter),
            # with respect to the posterior (ie measure information loss in the particle estimate,
            #   a sample from the posterior and a random sample, relative to the posterior)
            PosteriorEntropy = Float64[], 
            KLD              = Float64[], 
            KLD0             = Float64[], 
            KLDI             = Float64[],

            ## posterior density summary stats ##

            # posterior probability that predator is closer than 25, 50 and 100um
            PR40 = Float64[], 
            PR45 = Float64[], 
            PR50 = Float64[], 

            # quantiles of posterior density of distance to predator
            # 1%, 5%, 25% and 50% (median) (proximal side only, ie we care about how close the predator
            # might be, not how far away it might be)
            QP01 = Float64[], 
            QP05 = Float64[], 
            QP25 = Float64[], 
            QP50 = Float64[],

            # quantiles of angular distribution
            Qψ01 = Float64[], 
            Qψ05 = Float64[], 
            Qψ25 = Float64[], 
            Qψ50 = Float64[], 
            Qψ75 = Float64[], 
            Qψ95 = Float64[], 
            Qψ99 = Float64[],

            # M-cell Bayesian posterior belief of predator in patch
            MP = Float64[],

            ## Particle summary stats  ##

            # proportion of particles within 25, 50 and 100um
            NR40 = Float64[], 
            NR45 = Float64[], 
            NR50 = Float64[], 

            # quantiles of particle proximity, 1%, 5%, 25% and 50%
            # e.g. QN05 is range including closest 5% of particles
            QN01 = Float64[], 
            QN05 = Float64[], 
            QN25 = Float64[], 
            QN50 = Float64[], 

            # quantiles of particle direction error (from true heading to predator)
            # giving credibility intervals for direction
            # e.g. QΘ05 is left/anticlockwise limit of 5% credibility interval for direction to predator,
            #      and Q095 is right/clockwise limit
            QΘ01 = Float64[], 
            QΘ05 = Float64[], 
            QΘ25 = Float64[], 
            QΘ50 = Float64[], 
            QΘ75 = Float64[], 
            QΘ95 = Float64[], 
            QΘ99 = Float64[],

            # M-cell particle posterior belief of predator in patch
            MN = Float64[],
            
 ))

end # if DO_LOG_DATA


# SIMULATION LOOP STARTS HERE
tick()



for rep = 1:N_REPS

    # print trial parameters in terminal
    println("Replicate = ",  rep, 
            ", Likelihood particles = ", n_likelihood_particles, 
            ", Posterior particles = ", n_posterior_particles, 
            ", Prior bias  = ", posteriorDeathRate
            )

    # construct placozoans
    # HINT for debugging: 'prey' and 'predator' are declared global so that if
    #       you comment out 'placozoanStalker()' function (and its matching 'end')
    #       then they exist in the global REPL workspace when the program is interrupted by Ctrl-C.
    #       (otherwise they remain local within 'placozoanStalker()'  )
    global
    prey = Placozoan(prey_radius, prey_margin, prey_fieldrange,
        Nreceptors, sizeof_receptor, mat_radius,
        Ncrystals, sizeof_crystal, mat_radius,
        n_likelihood_particles, n_posterior_particles,
        posteriorDeathRate,  posterior_diffusion_coefficient, particle_collisionRadius, nFrames)
    prey.ω[] = prey_spinrate

    # predator 
    θ = 2.0*π * rand() # Random initial heading (from above) 
    predator_initial_position = (prey_radius + predator_radius + Δinit)*Point2(cos(θ), sin(θ))
    global
    predator = Placozoan(predator_radius, predator_margin, predator_fieldrange,
                            predator_initial_position, predator_speed, predator_brownian,
                            RGBA(.25, 0.1, 0.1, .35),
                            RGBA(.45, 0.1, 0.1, 0.2),
                            RGB(.95, 0.1, 0.1) )

    placozoanFieldstrength!(predator)
    Ereceptor_RF(prey, predator)
    initializePosteriorPDF(prey)       #  prior density on grid, copied to initial posterior 
    initializePosteriorParticles(prey) # initial posterior particles are sample from prior

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
        # prey.observer.prior[:,:] = prey.observer.posterior[:,:]
        radialSmooth(prey.observer.prior, prey_radius:mat_radius)           
    end

        

    # initialize particle filter
    #     initialize_particles(prey) # draw initial sample from prior

    if SHOW_ANIMATION

        # if !PLOT_ARRAYS    # not plotting likelihoods or posterior
        #     scene, layout = layoutscene(resolution=(WorldSize, WorldSize))
        #     left_panel =    layout[1, 1] = LAxis(scene, title="Placozoan", 
        #                                         backgroundcolor=colour_background )

        #     hidespines!(left_panel)
        #     hidedecorations!(left_panel)
        # else
        scene = Figure(resolution = ( Int(cround(3 * .75 * WorldSize)), Int(cround(.75 * WorldSize + 40)) ), backgroundcolor = :black)
        left_panel = scene[1,1] = Axis(scene, title="2D Reaction-Diffusion (Placozoan Model)",
            titlecolor = title_color, backgroundcolor=:black )
        middle_panel = scene[1, 2] = Axis(scene, title="Bayesian Observer", titlecolor = title_color, backgroundcolor=colour_background)
        right_panel = scene[1, 3] = Axis(scene, title="Likelihood", titlecolor = title_color, backgroundcolor=colour_background)





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
    # end


        # mat is a dark green disc in left panel (always present)
        mat_plt = poly!(left_panel,
            decompose(Point2f, Circle(Point2f(0., 0.), mat_radius)),
            color=colour_mat, strokewidth=0, strokecolor=:black)
        display(scene)

        mat_middle_plt = poly!(middle_panel,
            decompose(Point2f, Circle(Point2f(mat_radius, mat_radius), mat_radius)),
            color=RGBA(0.0, 0.0, 0.0, 0.0),
            strokewidth=0.5,
            strokecolor=RGB(0.75, 0.75, 0.75),
        )

        mat_right_plt = poly!(right_panel,
            decompose(Point2f, Circle(Point2f(mat_radius, mat_radius), mat_radius)),
            color=RGBA(0.0, 0.0, 0.0, 0.0),
            strokewidth=0.5,
            strokecolor=RGB(0.75, 0.75, 0.75),
        )

        # time and Δ report on plot
        text!(left_panel, @lift(string("T = ", Int64($t[]))), color = :white, 
                            position = (40-mat_radius,60-mat_radius), textsize = 16)
        # text!(left_panel, @lift(string("Δ = ", Int64(round(prey.observer.range[Int64($t[])]-prey.radius-predator.radius)))), color = :white, 
        #                   position = (40-mat_radius,35-mat_radius), textsize = 16)
        text!(left_panel, @lift(string("Δ = ", Int64(cround(prey.observer.Δ[Int64($t[])] ) ))), color = :white, 
                            position = (40-mat_radius,35-mat_radius), textsize = 16)                #prey.observer.range[i]

    # predator drawn using lift(..., node)
    # (predatorLocation does not depend explicitly on t, but this causes
    #  the plot to be updated when the node t changes)
    predator_plt = poly!(left_panel,
        lift(s -> decompose(Point2f, Circle(Point2f(predator.position[][1], predator.position[][2]),
            predator.radius)), t),
        color=predator.color, strokecolor=predator.edgecolor, strokewidth=.5)

        # plot likelihood particles (samples from likelihood)
        Lparticle_plt = scatter!(left_panel,prey.observer.Lparticle[:], 
            # prey.observer.Lparticle[1:prey.observer.nLparticles[], 1],
            # prey.observer.Lparticle[1:prey.observer.nLparticles[], 2],
            color=colour_likelihood, markersize=size_likelihood, strokewidth=0.1)

        # plot posterior particles
        Pparticle_plt = scatter!(left_panel, prey.observer.Pparticle[:],
            # prey.observer.Pparticle[1:prey.observer.nPparticles[], 1],
            # prey.observer.Pparticle[1:prey.observer.nPparticles[], 2],
            color=colour_posterior, markersize=size_posterior, strokewidth=0.1)

    # plot projection of likelihood particles into prey margin
    # nb this is a dummy plot
    # the correct particle locations are inserted before first plot
    observation_plt = scatter!(left_panel, prey.observer.Sparticle[:],
        # zeros(prey.observer.nLparticles[]),
        # zeros(prey.observer.nLparticles[]),
        color=:yellow, strokewidth=0, markersize=size_observation )

    # plot projection of posterior particles into prey margin
    # nb this is a dummy plot
    # the correct particle locations are inserted before first plot
    belief_plt = scatter!(left_panel, prey.observer.Bparticle[:], 
        # zeros(prey.observer.nPparticles[]),
        # zeros(prey.observer.nPparticles[]),
                            color=colour_posterior, strokewidth=0, markersize=size_belief)

    Likely_plt = plot!( right_panel,
        OffsetArrays.no_offset_view(prey.observer.likelihood), colorrange = (0.0, 1.25), colormap = :copper)  # :turku

    Posty_plt =  surface!( middle_panel , 1:WorldSize, 1:WorldSize,
        OffsetArrays.no_offset_view(prey.observer.posterior),  colormap = :magma)

    # PostContour_plt = contour!(right_panel, 1:WorldSize, 1:WorldSize,
    #     lift(u->u, Posty_plt[3]), levels = [1.0e-6, 1.0e-5, 1.0e-4], color = RGB(.35,.35,.35))

    predator_right_plt = poly!(right_panel,
        lift(s -> decompose(Point2f, Circle(Point2f(
            mat_radius + predator.position[][1], mat_radius + predator.position[][2]),
            predator.radius)), t ),
            color=RGBA(0.0, 0.0, 0.0, 0.0),
            strokecolor=predator.edgecolor,
            strokewidth=1.0,
        )

        predator_middle_plt = poly!( middle_panel,
            lift(s -> decompose(Point2f,
                Circle(Point2f(mat_radius + predator.position[][1], mat_radius + predator.position[][2]),
                predator.radius) ), t),
                color=RGBA(0.0, 0.0, 0.0, 0.0),
                strokecolor=predator.edgecolor,
                strokewidth=1.0    
            )

        prey_Lcopy_plt = poly!(middle_panel,
            decompose(Point2f, Circle(Point2f(mat_radius, mat_radius), prey.radius)),
            color=RGBA(0.0, 0.0, 0.0, 0.0), strokewidth=1.5, strokecolor=prey.gutcolor)

        preygut_Lcopy_plt = poly!(middle_panel,
            decompose(Point2f, Circle(Point2f(mat_radius, mat_radius), prey.radius - prey.marginwidth + 1)),
            color=prey.gutcolor, strokewidth=0.5, strokecolor=prey.gutcolor*.75)

        prey_Pcopy_plt = poly!(right_panel,
            decompose(Point2f, Circle(Point2f(mat_radius, mat_radius), prey.radius)),
            color=RGBA(0.0, 0.0, 0.0, 0.0), strokewidth=1.5, strokecolor=prey.gutcolor)

        preygut_Pcopy_plt = poly!(right_panel,
            decompose(Point2f, Circle(Point2f(mat_radius, mat_radius), prey.radius - prey.marginwidth + 1)),
            color=prey.gutcolor, strokewidth=0.5, strokecolor=prey.gutcolor*.75)



        # Prey
        prey_plt = poly!(left_panel,
            decompose(Point2f, Circle(Point2f(0., 0.), prey.radius)),
            color=prey.color, strokewidth=0.25, strokecolor=prey.gutcolor)

        preyGut_plt = poly!(left_panel,
            decompose(Point2f, Circle(Point2f(0., 0.), prey.gutradius)),
            color=prey.gutcolor, strokewidth=0.5, strokecolor=prey.gutcolor*.75)


        receptor_plt = scatter!(left_panel, prey.receptor.position[:], #prey.receptor.x, prey.receptor.y ,
            markersize=prey.receptor.size,
            color=[prey.receptor.openColor for i in 1:prey.receptor.N],
            strokecolor=:black, strokewidth=0.25)

        L_receptor_plt = scatter!(middle_panel, 
            (mat_radius + 1.0)*fill(Point2(1.0, 1.0), length(prey.receptor.position)) + prey.receptor.position,
            # mat_radius .+1 .+prey.receptor.x, mat_radius .+1 .+prey.receptor.y ,
            markersize=prey.receptor.size,
            color=[prey.receptor.openColor for i in 1:prey.receptor.N],
            strokecolor=:black, strokewidth=0.25)

        R_receptor_plt = scatter!(right_panel,  
            (mat_radius + 1.0)*fill(Point2(1.0, 1.0), length(prey.receptor.position)) + prey.receptor.position,
            #   mat_radius .+1 .+prey.receptor.x, mat_radius .+1 .+prey.receptor.y ,
            markersize=prey.receptor.size,
            color=[prey.receptor.openColor for i in 1:prey.receptor.N],
            strokecolor=:black, strokewidth=0.25)
            
        if PHOTORECEPTION
        crystal_plt = scatter!(left_panel, prey.photoreceptor.position[:], #prey.photoreceptor.x, prey.photoreceptor.y,
            markersize=prey.photoreceptor.size, marker=:diamond,
                color=[prey.photoreceptor.lightColor for i in 1:prey.photoreceptor.N],
                strokecolor=:black, strokewidth=0.25)
        end


        # reset axis limits (have been auto-adjusted by MakieLayout)
        xlims!(left_panel, -mat_radius, mat_radius)
        ylims!(left_panel, -mat_radius, mat_radius)

        xlims!(middle_panel, 0, WorldSize)
        ylims!(middle_panel, 0, WorldSize)
        xlims!(right_panel, 0, WorldSize)
        ylims!(right_panel, 0, WorldSize)


    end #SHOW_ANIMATION

    videoFileName = DataFileName * "_" * string(rep) * ".mkv"

    # VIDEO RECORDING ON/OFF
    # comment out ONE of the following 2 lines to generate video file or not
    # NB animation can be displayed while simulating, by setting SHOW_ANIMATION = true, without saving as video file,
    # SHOW_ANIMATION must be true for video to be recorded.
    #record(scene, videoFileName , framerate=30, compression = 0, 1:nFrames) do i     # simulate and write video file
    for i in 1:nFrames                                         # simulate without writing video file


        # predator takes a stochastic step toward prey
        # stalk() returns distance between edges at end of step
        prey.observer.Δ[i] = stalk(predator, prey, min_Δ)

        # electroreception
        if ELECTRORECEPTION
            electroreception(prey, predator)
        end  # electroreception

        # photoreception
        if PHOTORECEPTION
            photoreception(prey, predator)
            if SHOW_ANIMATION
                # set color of each receptor, indicating open or closed state
                crystalColor = [prey.photoreceptor.lightColor  for j in 1:prey.photoreceptor.N]
                crystalColor[findall(x -> x == 1, prey.photoreceptor.state)] .= prey.photoreceptor.darkColor
                crystal_plt.color[] = crystalColor
            end # SHOW_ANIMATION
        end # photoreception


        # inference
        # ELECTRORECEPTION & PHOTORECEPTION are bools 
        #   specifying whether the respective receptor states are included in the inference

        likelihood(prey, ELECTRORECEPTION, PHOTORECEPTION)  

        sample_likelihood(prey)     # random sample from likelihood

        bayesParticleUpdate(prey)   # Bayesian particle filter

        bayesArrayUpdate(prey)      # numerical sequential Bayes (benchmark)

        if SHOW_ANIMATION

            receptorColor = [prey.receptor.closedColor for j = 1:prey.receptor.N]
            receptorColor[findall(x -> x == 1, prey.receptor.state)] .= prey.receptor.openColor
            receptor_plt.color[] = L_receptor_plt.color[] = R_receptor_plt.color[] = receptorColor

            # update likelihood particle plot
            Lparticle_plt[1] = prey.observer.Lparticle
            #Lparticle_plt[1] = prey.observer.Lparticle[:]

            # update posterior particle plot
            Pparticle_plt[1] = prey.observer.Pparticle
            #Pparticle_plt[1] = prey.observer.Pparticle[:]

            # update observation particle plot
            observation_plt[1] = prey.observer.Sparticle

            # update observation particle plot
            belief_plt[1] = prey.observer.Bparticle


            Likely_plt[1] = mask .* OffsetArrays.no_offset_view(prey.observer.likelihood)
            Posty_plt[3] = mask .* OffsetArrays.no_offset_view(prey.observer.posterior)


        else
            # show progress in terminal
            # (so you know it's running even if SHOW_ANIMATION is false)
            print(".")
        end # SHOW_ANIMATION

        # record posterior entropy 
        prey.observer.PosteriorEntropy[i] = entropy(prey.observer.posterior)
        
        # # record distance between edges of predator and prey
        # getRange(prey, predator, i)

        # Kullback-Liebler divergence (information-distance between particle distibution and 'true' posterior)
        KLD!(prey, i)

        t[] = i  # update observable, causes redraw 
        notify(t)

        if LOG_DATA

            # MAP predator location            
            iMAP =  MAPlocation(prey)
            #scatter!(iMAP.+400, color = :blue)
            # println(iMAP[1], ", ", iMAP[2])

            # summary stats of particle distribution
            # NR = proportion of particles (estimated probabilty) that predator is closer than 25,50 & 100um
            # QN = quantiles of particle range, [0.005 0.025 0.25 0.5 0.75 0.975 0.995]
            #      giving 1%, 5% and 50% credibility intervals + median estimate
            # QΘ = quantiles of particle angle deviation from heading to predator (as above)
            # MN = M-cell's posterior belief that there is a predator in its patch 
            (NR, QN, Qθ, MN) = particleStats(prey, predator) 

            # corresponding stats for Bayesian observer
            (PR, QP, Qψ, MP) = observerStats(prey, predator) 

            # debug
            # println(QD, ", ", Dmin, ", ", Qθ, ", ", θmin, ", ", θmax)
            # sleep(2)

            # save data (see file open command for more detailed description of variables saved)
            CSV.write(DataFileName * ".csv",
                DataFrame(hcat(rep, 
                
                # trial parameters
                n_likelihood_particles, 
                n_posterior_particles,  posteriorDeathRate,
                                        # range and location 
                prey.observer.Δ[i], predator.position[][1], predator.position[][2], iMAP[1], iMAP[2],

                # entropy/information in particle filter
                prey.observer.PosteriorEntropy[i], prey.observer.KLD[i], prey.observer.KLD0[i], prey.observer.KLDI[i], 
                
                # summary stats of posterior probability 
                PR..., QP, Qψ, MP, 

                # summary stats of particle distribution
                NR..., QN..., Qθ..., MN), 
                
                :auto), header=false, append=true)

            end # if LOG_DATA   

        end # frame

        # print timing information to terminal
        println()
        laptimer()

        t = 1

    end # for rep

end # placozoanStalker function


# run it
placozoanStalker()





   
