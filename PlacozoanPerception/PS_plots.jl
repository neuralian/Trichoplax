using CairoMakie, ColorSchemes, Colors
include("makieTheme5.jl") 

#############################################################
# FIGURE 1: Entropy and K-L Divergence 
fig1 = Figure(resolution = (800,600))
fig1[1,1] = ax = Axis(fig1, xticks = [25 50 75 100 125], 
    title = "Particle Filter Performance",
    xlabel = "predator range /μm",
    ylabel = "K-L Divergence (bits)")


lines!(Rgrid,KLDS[:,iQ50], color = "#56B4E9")
# lines!(ax,KLDS[:,iQ25], color = :red)
# lines!(ax,KLDS[:,iQ75], color = :red)
band!(Rgrid,KLDS[:,iQ25] ,KLDS[:,iQ75], color = ("#56B4E9", 0.25) )

lines!(Rgrid,KLD0S[:,iQ50], color = "#009E73")
band!(Rgrid,KLD0S[:,iQ25] ,KLD0S[:,iQ75], color = ("#009E73", 0.25) )

lines!(Rgrid,KLDIS[:,iQ50], color = "#E69F00")
band!(Rgrid,KLDIS[:,iQ25] ,KLDIS[:,iQ75], color = ("#E69F00", 0.25) )

lines!(Rgrid,-PENTS[:,iQ50], color = "#E84646")
band!(Rgrid,-PENTS[:,iQ75] ,-PENTS[:,iQ25], color = ("#E84646", 0.25) )

ax.xreversed = true
#xlims!(25, 125)
display(fig1);


save("KLD_fig.png", fig1, px_per_unit = 3 )

###############################################################
# FIGURE 2: Probability of Predator in range
fig2 = Figure(resolution = (800,600))
fig2[1,1] = ax = Axis(fig2, xticks = [25 50 100], 
    yticks = [0 0.25 0.5 0.75 1.0],
    title = "Probability of predator within 50μm",
    xlabel = "predator range /μm",
    ylabel = "probability")

#lines!(Rgrid, PR25S[:, iQ50])
S = sum(PR50S[:, iQ50])/100.
lines!(Rgrid, PR50S[:, iQ50]./S, color = "#E69F00")
band!(Rgrid, PR50S[:, iQ25]./S, PR50S[:, iQ75]./S, color = ("#E69F00", .1) )
#lines!(Rgrid, PR100S[:, iQ50])

# lines!(Rgrid, NR25S[:, iQ50], color = :blue)
lines!(Rgrid, NR50S[:, iQ50], color = "#56B4E9" )
band!(Rgrid, NR50S[:, iQ25], NR50S[:, iQ75], color = ("#56B4E9", .1))
# lines!(Rgrid, NR100S[:, iQ50], color = :blue)


ax.xreversed = true

display(fig2);


save("PR_fig.png", fig2, px_per_unit = 3 )

###############################################################
# FIGURE 3: Quantiles of range Distribution
fig3 = Figure(resolution = (800,600))
ax = Axis(fig3, #xticks = [25 50 100], 
   # yticks = [0 0.25 0.5 0.75 1.0],
    title = "Median (.05,.50) Credibility Interval for Distance to Predator",
    xlabel = "True Range /μm",
    ylabel = "Inferred Range")


BayesLine = lines!(Rgrid, QP50S[:,iQ50], color = "#E69F00", linewidth = 0.5)
band!(Rgrid, QP05S[:,iQ50], QP50S[:,iQ50], color = ("#E69F00", .1) )
lines!(Rgrid, QP05S[:,iQ50], color = "#E69F00", linewidth = 0.5)

ParticleLine = lines!(Rgrid, QN50S[:,iQ50], color = "#56B4E9", linewidth = 0.5)
band!(Rgrid, QN05S[:,iQ50], QN50S[:,iQ50], color = ("#56B4E9", .1) )
lines!(Rgrid, QN05S[:,iQ50], color = "#56B4E9", linewidth = 0.5)

lines!(Rgrid, Rgrid, color = "#E84646" )

leg3 = Legend(fig3, [BayesLine ParticleLine], 
                    [" Numerical Sequential Bayes", "Placozoan Particle Filter"], 
            halign = :right, valign = :top, markersize = 4)

fig3[1,1] = ax
fig3[1,1] = leg3

ax.xreversed = true

display(fig3);


save("QP_fig.png", fig3, px_per_unit = 3 )

###############################################################
# FIGURE 4: Quantiles of angle distribution
fig4 = Figure(resolution = (800,600))
ax = Axis(fig4, #xticks = [25 50 100], 
   # yticks = [0 0.25 0.5 0.75 1.0],
    title = "Credibility Interval for Heading to Predator",
    xlabel = "True Range /μm",
    ylabel = "Inferred Angle")

lines!(Rgrid, QΘ50S[:, iQ50], color = "#E69F00", linewidth = 1.0)
band!(Rgrid, QΘ25S[:, iQ50], QΘ75S[:, iQ50], color = ("#E69F00", 0.1) )

fig4[1,1] = ax
ax.xreversed = true

save("Qθ_fig.png", fig4, px_per_unit = 3 )

display(fig4);