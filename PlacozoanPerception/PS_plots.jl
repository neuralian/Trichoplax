# Run PS_readDataFrame before this

using CairoMakie, ColorSchemes, Colors
include("makieThemeBlue.jl") 


BayesColor = "#E69F00"
ParticleColor = "#56B4E9"
reddishColor =  "#E84646"
greenishColor = "#009E73"
thickLineW = 3
lineW = 1.5
traceW = 0.5
bandAlpha = 0.5

(Nrows, Ncols) = size(KLD)
ieg = rand(1:Ncols)  # pick a random trial to show
Neg = 8  # number of random trials to show


if 1==2
#############################################################
# FIGURE 1: Entropy
fig1 = Figure(resolution = (800,600))
ax1 = Axis(fig1, xticks = [25 50 75 100 125], 
    title = "Information/Uncertainty about Location",
    xlabel = "Separation (μm)",
    ylabel = "Entropy / Information (bits)")


# # plaocozoan posterior entropy
# pLine = lines!(Rgrid,-KLDS[:,iQ50], color =ParticleColor)
# [ lines!(Rgrid, -KLDS[:, rand(1:Ncols)], color = ParticleColor, linewidth = traceW) for _ in 1:Neg]
# band!(Rgrid,-KLDS[:,iQ75] ,-KLDS[:,iQ25], color = (ParticleColor, bandAlpha) )

# K-L divergence of uniform random from Bayesian posterior 
zLine = lines!(Rgrid,-KLD0S[:,iQ50], color = greenishColor, linewidth = thickLineW)
band!(Rgrid,-KLD0S[:,iQ75] ,-KLD0S[:,iQ25], color = (greenishColor, bandAlpha) )
[ lines!(Rgrid, -KLD0[:, rand(1:Ncols)], color = greenishColor, linewidth = traceW) for _ in 1:Neg]


# iLine = lines!(Rgrid,-KLDIS[:,iQ50], color = BayesColor)
# band!(Rgrid,-KLDIS[:,iQ75] ,-KLDIS[:,iQ25], color = (BayesColor, bandAlpha) )

# Entropy of Bayesian Posterior
eLine = lines!(Rgrid,PENTS[:,iQ50], color = BayesColor, linewidth = thickLineW)
band!(Rgrid,PENTS[:,iQ75] ,PENTS[:,iQ25], color = (BayesColor, bandAlpha) )
[ lines!(Rgrid, PENT[:, rand(1:Ncols)], color = BayesColor, linewidth = traceW) for _ in 1:Neg]

ax1.xreversed = true

leg1 = Legend(fig1, [eLine zLine ], 
[" Posterior Entropy",  " Uniform Divergence"], 
            halign = :right, valign = :bottom, markersize = 4, labelsize = 12)

fig1[1,1] = ax1
fig1[1,1] = leg1


#xlims!(25, 125)
display(fig1);

save("Entropy_fig.png", fig1, px_per_unit = 3 )


#############################################################
# FIGURE 1a: K-L Divergence 
fig1a = Figure(resolution = (800,600))
ax1a = Axis(fig1a, xticks = [25 50 75 100 125], 
    title = "Information Loss (K-L Divergence)",
    xlabel = "Separation (μm)",
    ylabel = "Entropy / Information (bits)")


# plaocozoan posterior entropy
pLine = lines!(Rgrid,-KLDS[:,iQ50], color =ParticleColor, linewidth = thickLineW)
[ lines!(Rgrid, -KLD[:, rand(1:Ncols)], color = ParticleColor, linewidth = traceW) for _ in 1:Neg]
band!(Rgrid,-KLDS[:,iQ75] ,-KLDS[:,iQ25], color = (ParticleColor, bandAlpha) )

# # K-L divergence of uniform random from Bayesian posterior 
# zLine = lines!(Rgrid,-KLD0S[:,iQ50], color = greenishColor)
# band!(Rgrid,-KLD0S[:,iQ75] ,-KLD0S[:,iQ25], color = (greenishColor, bandAlpha) )
# [ lines!(Rgrid, -KLD0S[:, rand(1:Ncols)], color = greenishColor, linewidth = traceW) for _ in 1:Neg]

# ideal particle filter
iLine = lines!(Rgrid,-KLDIS[:,iQ50], color = BayesColor)
[ lines!(Rgrid, -KLDI[:, rand(1:Ncols)], color = BayesColor, linewidth = traceW) for _ in 1:Neg]
band!(Rgrid,-KLDIS[:,iQ75] ,-KLDIS[:,iQ25], color = (BayesColor, bandAlpha) )

# eLine = lines!(Rgrid,PENTS[:,iQ50], color = reddishColor)
# band!(Rgrid,PENTS[:,iQ75] ,PENTS[:,iQ25], color = (reddishColor, bandAlpha) )

ax1a.xreversed = true

leg1a = Legend(fig1a, [pLine iLine], 
[" Placozoan", " Bayesian Particle Filter"], 
            halign = :left, valign = :bottom, markersize = 4, labelsize = 12)

fig1[1,1] = ax1a
fig1[1,1] = leg1a

save("KLD_fig.png", fig1a, px_per_unit = 3 )

#xlims!(25, 125)
display(fig1a);







###############################################################
# FIGURE 2: Probability of Predator in range
fig2 = Figure(resolution = (800,600))
ax2 = Axis(fig2,
    yticks = [0 0.25 0.5 0.75 1.0],
    title = "Proximity Detection",
    xlabel = "Separation (μm)",
    ylabel = "Posterior Probability of Separation <50μm",
    xticks = [25 50 75 100 125])

#lines!(Rgrid, PR25S[:, iQ50])
# S = sum(PR50S[:, iQ50])/100.
band!(Rgrid, PR50S[:, iQ25], PR50S[:, iQ75], color = (BayesColor, bandAlpha) )
BayesLine = lines!(Rgrid, PR50S[:, iQ50], color = BayesColor, linewidth = thickLineW)
[ lines!(Rgrid, PR50[:, rand(1:Ncols)], color = BayesColor, linewidth = traceW) for _ in 1:Neg]

# lines!(Rgrid, NR25S[:, iQ50], color = :blue)
# S1  = sum(NR50S[:, iQ50])/100.
band!(Rgrid, NR50S[:, iQ25], NR50S[:, iQ75], color = (ParticleColor, bandAlpha))
ParticleLine = lines!(Rgrid, NR50S[:, iQ50], color = ParticleColor, linewidth = thickLineW)
[ lines!(Rgrid, NR50[:, rand(1:Ncols)], color = ParticleColor, linewidth = traceW) for _ in 1:Neg]

leg2 = Legend(fig2, [BayesLine ParticleLine], 
                    [" Bayesian", " Placozoan"], 
            halign = :left, valign = :top, markersize = 4, labelsize = 12)

fig2[1,1] = ax2
fig2[1,1] = leg2

ax2.xreversed = true

display(fig2);


save("PR_fig.png", fig2, px_per_unit = 3 )

end


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
BayesLine = lines!(Rgrid, QP50S[:,iQ50], color = BayesColor, linewidth = thickLineW)
[ lines!(Rgrid, QP50[:, rand(1:Ncols)], color = BayesColor, linewidth = traceW) for _ in 1:Neg]


band!(Rgrid, QN05S[:,iQ50], QN50S[:,iQ50], color = (ParticleColor, bandAlpha) )
ParticleLine = lines!(Rgrid, QN50S[:,iQ50], color = ParticleColor, linewidth = thickLineW)
[ lines!(Rgrid, QN50[:, rand(1:Ncols)], color = ParticleColor, linewidth = traceW) for _ in 1:Neg]

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
    xlabel = "Separation (μm)",
    ylabel = "Error in Inferred Bearing (degrees)",
    yticks = [-45 -5 5 45],
    xticks = [25 50 75 100 125])



band!(Rgrid, QΘ25S[:, iQ50], QΘ75S[:, iQ50], color = (ParticleColor, bandAlpha) )
ParticleLine4 = [ lines!(Rgrid, QΘ50[:, rand(1:Ncols)], color = ParticleColor, linewidth = traceW) for _ in 1:Neg]

band!(Rgrid, Qψ25S[:, iQ50], Qψ75S[:, iQ50], color = (BayesColor,bandAlpha) )
BayesLine4=[ lines!(Rgrid, Qψ50[:, rand(1:Ncols)], color = BayesColor, linewidth = traceW) for _ in 1:Neg]

leg4 = Legend(fig4, [BayesLine4[1] ParticleLine4[1]], 
                    [" Bayesian", " Placozoan"], 
            halign = :right, valign = :top, markersize = 4, labelsize = 12)

fig4[1,1] = ax4
fig4[1,1] = leg4
ax4.xreversed = true

save("Qθ_fig.png", fig4, px_per_unit = 3 )

display(fig4);




###############################################################
# FIGURE 5: Mauthner Cell
fig5 = Figure(resolution = (800,600))
ax5 = Axis(fig5, #xticks = [25 50 100], 
   # yticks = [0 0.25 0.5 0.75 1.0],
    title = "M-cell",
    xlabel = "Separation (μm)",
    ylabel = "Probability",
    xticks = [25 50 75 100 125])


SC = 1.0
band!(Rgrid, MCPS[:, iQ25], MCPS[:, iQ75], color = (BayesColor, bandAlpha) )
[ lines!(Rgrid, MCP[:, rand(1:Ncols)], color = BayesColor, linewidth = traceW) for _ in 1:Neg]
BayesLine5=lines!(Rgrid, MCPS[:, iQ50], color = BayesColor, linewidth = thickLineW)

band!(Rgrid, SC*MCNS[:, iQ25], SC*MCNS[:, iQ75], color = (ParticleColor, bandAlpha) )
[ lines!(Rgrid, MCN[:, rand(1:Ncols)], color = ParticleColor, linewidth = traceW) for _ in 1:Neg]
lines!(Rgrid, SC*MCNS[:, iQ50], color = ParticleColor, linewidth = lineW)



leg5 = Legend(fig5, [BayesLine5 ParticleLine5], 
                    [" Bayesian", "Placozoan"], 
            halign = :left, valign = :top, markersize = 4, labelsize = 12)

fig5[1,1] = ax5
fig5[1,1] = leg5
ax5.xreversed = true

save("MC_fig.png", fig5, px_per_unit = 3 )

display(fig5);

