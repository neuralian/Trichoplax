using CairoMakie, ColorSchemes, Colors
include("makieThemeBlue.jl") 

BayesColor = "#E69F00"
ParticleColor = "#56B4E9"
reddishColor =  "#E84646"
greenishColor = "#009E73"
lineW = 1.5
bandAlpha = 0.35

(Nrows, Ncols) = size(KLD)
ieg = rand(1:Ncols)  # pick a random trial to show


#############################################################
# FIGURE 1: Entropy and K-L Divergence 
fig1 = Figure(resolution = (800,600))
ax1 = Axis(fig1, xticks = [25 50 75 100 125], 
    title = "Location Uncertainty",
    xlabel = "Separation (μm)",
    ylabel = "Entropy / Information (bits)")


pLine = lines!(Rgrid,-KLDS[:,iQ50], color =ParticleColor)
# lines!(ax,KLDS[:,iQ25], color = :red)
# lines!(ax,KLDS[:,iQ75], color = :red)
band!(Rgrid,-KLDS[:,iQ75] ,-KLDS[:,iQ25], color = (ParticleColor, bandAlpha) )

zLine = lines!(Rgrid,-KLD0S[:,iQ50], color = greenishColor)
band!(Rgrid,-KLD0S[:,iQ75] ,-KLD0S[:,iQ25], color = (greenishColor, bandAlpha) )

iLine = lines!(Rgrid,-KLDIS[:,iQ50], color = BayesColor)
band!(Rgrid,-KLDIS[:,iQ75] ,-KLDIS[:,iQ25], color = (BayesColor, bandAlpha) )

eLine = lines!(Rgrid,PENTS[:,iQ50], color = reddishColor)
band!(Rgrid,PENTS[:,iQ75] ,PENTS[:,iQ25], color = (reddishColor, bandAlpha) )

ax1.xreversed = true

leg1 = Legend(fig1, [eLine pLine zLine iLine], 
[" Posterior Entropy", " Placozoan Divergence ", " Uniform Divergence", " Optimal Divergence"], 
            halign = :left, valign = :bottom, markersize = 4, labelsize = 12)

fig1[1,1] = ax1
fig1[1,1] = leg1


#xlims!(25, 125)
display(fig1);


save("KLD_fig.png", fig1, px_per_unit = 3 )



###############################################################
# FIGURE 2: Probability of Predator in range
fig2 = Figure(resolution = (800,600))
ax2 = Axis(fig2,
    yticks = [0 0.25 0.5 0.75 1.0],
    title = "Proximity Detection",
    xlabel = "True Separation (μm)",
    ylabel = "Posterior Probability of Separation <50μm",
    xticks = [25 50 75 100 125])

#lines!(Rgrid, PR25S[:, iQ50])
# S = sum(PR50S[:, iQ50])/100.
# S=1
band!(Rgrid, PR50S[:, iQ25]./S, PR50S[:, iQ75]./S, color = (BayesColor, bandAlpha) )
lines!(Rgrid, PR50S[:, iQ50]./S, color = BayesColor, linewidth = lineW)
#lines!(Rgrid, PR100S[:, iQ50])

# lines!(Rgrid, NR25S[:, iQ50], color = :blue)
# S1  = sum(NR50S[:, iQ50])/100.
# S1 = 1
band!(Rgrid, NR50S[:, iQ25]/S1, NR50S[:, iQ75]/S1, color = (ParticleColor, bandAlpha))
lines!(Rgrid, NR50S[:, iQ50]/S1, color = ParticleColor, linewidth = lineW)
# lines!(Rgrid, NR100S[:, iQ50], color = :blue)

leg2 = Legend(fig2, [BayesLine ParticleLine], 
                    [" Bayesian", " Placozoan"], 
            halign = :left, valign = :top, markersize = 4, labelsize = 12)

fig2[1,1] = ax2
fig2[1,1] = leg2

ax2.xreversed = true

display(fig2);


save("PR_fig.png", fig2, px_per_unit = 3 )


###############################################################
# FIGURE 3: Quantiles of range Distribution
fig3 = Figure(resolution = (800,600))
ax3 = Axis(fig3, #xticks = [25 50 100], 
   # yticks = [0 0.25 0.5 0.75 1.0],
    title = "Separation Credibility",
    xlabel = "True Separation (μm)",
    ylabel = "Inferred Separation (μm)",
    xticks = [25 50 75 100 125],
    yticks = [25 50 75 100 125 150 175 200])



band!(Rgrid, QP05S[:,iQ50], QP50S[:,iQ50], color = (BayesColor, bandAlpha) )
#lines!(Rgrid, QP05S[:,iQ50], color = BayesColor, linewidth = lineW/2)
BayesLine = lines!(Rgrid, QP05[:,ieg], color = BayesColor, linewidth = lineW)


band!(Rgrid, QN05S[:,iQ50], QN50S[:,iQ50], color = (ParticleColor, bandAlpha) )
#lines!(Rgrid, QN05S[:,iQ50], color = ParticleColor, linewidth = lineW/2)
ParticleLine = lines!(Rgrid, QN50[:,ieg], color = ParticleColor, linewidth = lineW)

lines!(Rgrid, Rgrid, color = reddishColor, linewidth = lineW )

leg3 = Legend(fig3, [BayesLine ParticleLine], 
                    [" Bayesian", " Placozoan"], 
            halign = :right, valign = :top, markersize = 4, labelsize = 12)

fig3[1,1] = ax3
fig3[1,1] = leg3

ax3.xreversed = true

display(fig3);


save("QP_fig.png", fig3, px_per_unit = 3 )



###############################################################
# FIGURE 4: Quantiles of angle distribution
fig4 = Figure(resolution = (800,600))
ax4 = Axis(fig4, #xticks = [25 50 100], 
   # yticks = [0 0.25 0.5 0.75 1.0],
    title = "Direction Credibility",
    xlabel = "True Separation (μm)",
    ylabel = "Error in Inferred Bearing (degrees)",
    yticks = [-45 -5 5 45],
    xticks = [25 50 75 100 125])



band!(Rgrid, QΘ25S[:, iQ50], QΘ75S[:, iQ50], color = (ParticleColor, bandAlpha) )
ParticleLine4 = lines!(Rgrid, QΘ50[:, ieg], color = ParticleColor, linewidth = lineW)

band!(Rgrid, Qψ25S[:, iQ50], Qψ75S[:, iQ50], color = (BayesColor,bandAlpha) )
BayesLine4 = lines!(Rgrid, Qψ50[:, ieg], color = BayesColor, linewidth = lineW)

leg4 = Legend(fig4, [BayesLine4 ParticleLine4], 
                    [" Bayesian", " Placozoan"], 
            halign = :right, valign = :top, markersize = 4, labelsize = 12)

fig4[1,1] = ax4
fig4[1,1] = leg4
ax4.xreversed = true

save("Qθ_fig.png", fig4, px_per_unit = 3 )

display(fig4);




###############################################################
# FIGURE 5: Cellular proximity detector
fig5 = Figure(resolution = (800,600))
ax5 = Axis(fig5, #xticks = [25 50 100], 
   # yticks = [0 0.25 0.5 0.75 1.0],
    title = "Cellular Integration",
    xlabel = "Separation (μm)",
    ylabel = "Probability",
    xticks = [25 50 75 100 125])



band!(Rgrid, CPDS[:, iQ25], CPDS[:, iQ75], color = (ParticleColor, bandAlpha) )
BayesLine5 = lines!(Rgrid, CPD[:, ieg], color = greenishColor, linewidth = lineW)
lines!(Rgrid, CPDS[:, iQ50], color = ParticleColor, linewidth = lineW)


# leg5 = Legend(fig5, [BayesLine5 ], 
#                     [" Bayesian"], 
#             halign = :left, valign = :top, markersize = 4, labelsize = 12)

fig5[1,1] = ax5
fig5[1,1] = leg5
ax5.xreversed = true

save("CPD_fig.png", fig5, px_per_unit = 3 )

display(fig5);