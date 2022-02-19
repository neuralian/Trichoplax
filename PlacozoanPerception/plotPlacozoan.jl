# Run PS_readDataFrame before this

using CairoMakie, ColorSchemes, Colors
include("makieTheme1.jl") 


BayesBG = "#8793fa" # "#E69F00"
BayesFG = "#2d3bba"
ParticleBG = "#fa8787"  # #56B4E9"
ParticleFG = "#c92020"
reddishColor =  "#009E73"
greenishFG = "#009E73"
greenishBG = "#6df7d2"
thickLineW = 3
lineW = 1.5
traceW = 0.25
bandAlpha = 0.5

(Nrows, Ncols) = (Nframes, Nreps)
ieg = rand(1:Ncols)  # pick a random trial to show
Neg = 8  # number of random trials to show


Rgrid = Range
dticks = vec([5 25 50 75 100 125 150 175]) # xticks for distance

#############################################################
# FIGURE 1: Entropy
fig1 = Figure(resolution = (800,600))
ax1 = Axis(fig1[1,1], 
     xticks = dticks,
     yticks = vec([10 12 14 16 18 20]), 
    title = "Bayesian Uncertainty about Distance to Predator",
    xlabel = "Separation (μm)",
    ylabel = "Entropy (bits)")


# # plaocozoan posterior entropy
# pLine = lines!(Rgrid,-KLDS[:,iQ50], color =ParticleBG)
# [ lines!(Rgrid, -KLDS[:, rand(1:Ncols)], color = ParticleBG, linewidth = traceW) for _ in 1:Neg]
# band!(Rgrid,-KLDS[:,iQ75] ,-KLDS[:,iQ25], color = (ParticleBG, bandAlpha) )

# # K-L divergence of uniform random from Bayesian posterior 
# zLine = lines!(Rgrid,KLD0S[:,iQ50], color = greenishFG, linewidth = thickLineW)
# band!(Rgrid,KLD0S[:,iQ75] ,KLD0S[:,iQ25], color = (greenishBG, bandAlpha) )
# [ lines!(Rgrid, KLD0[:, rand(1:Ncols)], color = greenishFG, linewidth = traceW) for _ in 1:Neg]


# iLine = lines!(Rgrid,-KLDIS[:,iQ50], color = BayesBG)
# band!(Rgrid,-KLDIS[:,iQ75] ,-KLDIS[:,iQ25], color = (BayesBG, bandAlpha) )

# Entropy of Bayesian Posterior
eLine = lines!(Rgrid,PosteriorEntropy_sorted[:,iQ50], color = BayesFG, linewidth = thickLineW)
band!(Rgrid,PosteriorEntropy_sorted[:,iQ75] ,PosteriorEntropy_sorted[:,iQ25], color = (BayesBG, bandAlpha) )
[ lines!(Rgrid, PosteriorEntropy[:, rand(1:Ncols)], color = BayesFG, linewidth = traceW) for _ in 1:Neg]

ax1.xreversed = true

# leg1 = Legend(fig1, [eLine zLine ], 
# [" Bayesian Observer",  " Uniform Particle Distribution"], 
#             halign = :right, valign = :bottom, markersize = 4, labelsize = 12)

fig1[1,1] = ax1
#fig1[1,1] = leg1
ylims!(ax1,10, 20)

#xlims!(25, 125)
display(fig1);

save("Entropy_fig.png", fig1, px_per_unit = 3 )


#############################################################
# FIGURE 1a: K-L Divergence 
fig1a = Figure(resolution = (800,600))
ax1a = Axis(fig1a, 
    xticks = dticks, 
    yticks = vec([.25 .5 1.0 2 4 8 10]),
    title = "Information Loss relative to Bayesian Observer",
    xlabel = "Separation (μm)",
    ylabel = "Kullback-Liebler Divergence (bits)")




# # K-L divergence of uniform random from Bayesian posterior 
# K-L divergence of uniform random from Bayesian posterior 
zLine = lines!(Rgrid,KLD_random_sorted[:,iQ50], color = greenishFG, linewidth = thickLineW)
band!(Rgrid,KLD_random_sorted[:,iQ75] , KLD_random_sorted[:,iQ25], color = (greenishBG, bandAlpha) )
[ lines!(Rgrid, KLD_random[:, rand(1:Ncols)], color = greenishFG, linewidth = traceW) for _ in 1:Neg]

# ideal particle filter
iLine = lines!(Rgrid,KLD_sample_sorted[:,iQ50], color = BayesFG, linewidth = thickLineW)
[ lines!(Rgrid, KLD_sample[:, rand(1:Ncols)], color = BayesFG, linewidth = traceW) for _ in 1:Neg]
band!(Rgrid,KLD_sample_sorted[:,iQ75] ,KLD_sample_sorted[:,iQ25], color = (BayesBG, bandAlpha) )

# plaocozoan posterior entropy
pLine = lines!(Rgrid,KLD_particles_sorted[:,iQ50], color =ParticleFG, linewidth = thickLineW)
[ lines!(Rgrid, KLD_particles[:, rand(1:Ncols)], color = ParticleFG, linewidth = traceW) for _ in 1:Neg]
band!(Rgrid,KLD_particles_sorted[:,iQ75] , KLD_particles_sorted[:,iQ25], color = (ParticleBG, bandAlpha) )
ax1a.xreversed = true

leg1a = Legend(fig1a, [pLine iLine zLine], 
[" Placozoan Model", " Posterior Samples", " Uniform Samples"], 
            halign = :left, valign = :top, markersize = 4, labelsize = 12)

fig1[1,1] = ax1a
fig1[1,1] = leg1a

ylims!(ax1a, 0, 10)



#xlims!(25, 125)
display(fig1a);

save("KLD_fig.png", fig1a, px_per_unit = 3 )





###############################################################
# FIGURE 2: Probability of Predator in range
fig2 = Figure(resolution = (800,600))
ax2 = Axis(fig2,
    yticks = vec([0 0.25 0.5 0.75 1.0]),
    title = "Proximity Detection (Probability of Predator closer than 45μm)",
    xlabel = "Separation (μm)",
    ylabel = "Posterior Probability",
    xticks = dticks)

#lines!(Rgrid, PR25S[:, iQ50])
# S = sum(PR50S[:, iQ50])/100.
band!(Rgrid, P50_sorted[:, iQ25], P50_sorted[:, iQ75], color = (BayesBG, bandAlpha) )
BayesLine = lines!(Rgrid, P50_sorted[:, iQ50], color = BayesFG, linewidth = thickLineW)
[ lines!(Rgrid, P50[:, rand(1:Ncols)], color = BayesFG, linewidth = traceW) for _ in 1:Neg]

# lines!(Rgrid, NR25S[:, iQ50], color = :blue)
# S1  = sum(NR50S[:, iQ50])/100.
band!(Rgrid, N40_sorted[:, iQ25], N40_sorted[:, iQ75], color = (ParticleBG, bandAlpha))
ParticleLine = lines!(Rgrid, N40_sorted[:, iQ50], color = ParticleFG, linewidth = thickLineW)
[ lines!(Rgrid, N40[:, rand(1:Ncols)], color = ParticleFG, linewidth = traceW) for _ in 1:Neg]

leg2 = Legend(fig2, [BayesLine ParticleLine], 
                    [" Bayesian Observer", " Placozoan Model"], 
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
    title = "(0.05 - 0.5) Credibility Interval for Distance to Predator",
    xlabel = "True Separation (μm)",
    ylabel = "Inferred Separation (μm)",
    xticks = dticks,
    yticks = vec([0 25 50 75 100 125 150 175 200])  )



band!(Rgrid, QP05_sorted[:,iQ50], QP50_sorted[:,iQ50], color = (BayesBG, bandAlpha) )
BayesLine = lines!(Rgrid, QP50_sorted[:,iQ50], color = BayesFG, linewidth = thickLineW)
[ lines!(Rgrid, QP50[:, rand(1:Ncols)], color = BayesFG, linewidth = traceW) for _ in 1:Neg]


band!(Rgrid, QN05_sorted[:,iQ50], QN50_sorted[:,iQ50], color = (ParticleBG, bandAlpha) )
ParticleLine = lines!(Rgrid, QN50_sorted[:,iQ50], color = ParticleFG, linewidth = thickLineW)
[ lines!(Rgrid, QN50[:, rand(1:Ncols)], color = ParticleFG, linewidth = traceW) for _ in 1:Neg]

lines!(Rgrid, Rgrid, color = reddishColor, linewidth = thickLineW )

leg3 = Legend(fig3, [BayesLine ParticleLine], 
                    [" Bayesian Observer", " Placozoan Model"], 
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
    title = "Interquartile Range of Direction Error Credibility",
    xlabel = "Separation (μm)",
    ylabel = "Error in Inferred Bearing (degrees)",
    yticks = vec([-45 -5 5 45]),
    xticks = dticks)



band!(Rgrid, QΘ25_sorted[:, iQ50], QΘ75_sorted[:, iQ50], color = (ParticleBG, bandAlpha) )
ParticleLine4 = [ lines!(Rgrid, QΘ50[:, rand(1:Ncols)], color = ParticleFG, linewidth = traceW) for _ in 1:Neg]

band!(Rgrid, Qψ25_sorted[:, iQ50], Qψ75_sorted[:, iQ50], color = (BayesBG,bandAlpha) )
BayesLine4=[ lines!(Rgrid, Qψ50[:, rand(1:Ncols)], color = BayesFG, linewidth = traceW) for _ in 1:Neg]

leg4 = Legend(fig4, [BayesLine4[1] ParticleLine4[1]], 
                    [" Bayesian Observer", " Placozoan Model"], 
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
    title = "M-cell Activation",
    xlabel = "Separation (μm)",
    ylabel = "Probability",
    xticks = dticks)


SC = 1.0
band!(Rgrid, MP_sorted[:, iQ25], MP_sorted[:, iQ75], color = (BayesBG, bandAlpha) )
[ lines!(Rgrid, MP[:, rand(1:Ncols)], color = BayesFG, linewidth = traceW) for _ in 1:Neg]
BayesLine5=lines!(Rgrid, MP_sorted[:, iQ50], color = BayesFG, linewidth = thickLineW)

band!(Rgrid, SC*MN_sorted[:, iQ25], SC*MN_sorted[:, iQ75], color = (ParticleBG, bandAlpha) )
[ lines!(Rgrid, MN[:, rand(1:Ncols)], color = ParticleFG, linewidth = traceW) for _ in 1:Neg]
ParticleLine5=lines!(Rgrid, SC*MN_sorted[:, iQ50], color = ParticleFG, linewidth = thickLineW)



leg5 = Legend(fig5, [BayesLine5 ParticleLine5], 
                    [" Bayesian Observer", " Placozoan Model"], 
            halign = :left, valign = :top, markersize = 4, labelsize = 12)

fig5[1,1] = ax5
fig5[1,1] = leg5
ax5.xreversed = true

save("MC_fig.png", fig5, px_per_unit = 3 )

display(fig5);

