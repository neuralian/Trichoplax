# importPlacozoansData.jl
# multi-file version of importPlacozoanData.jl
# Reads .csv files created by PlacozoanStalker 
# Creates arrays convenient for plotPlacozoanData.jl
# and does some preprocessing/filtering
# MGP 2021-2022


using CSV
using DataFrames
using Statistics
using Dierckx        # spline interpolation


#include("makieTheme5.jl")

# read data file
fileNames = ["PlacozoanStalkerV2.0 2022-02-23 18.11.csv"
             ] 
RawData = nothing
Nreps = 0
Nframes = 0
# nb dropmissing allows deleting incomplete reps from crashed runs
# (deleted rows are imported as rows of missings)
for fileName in fileNames
    
    global RawData, Nreps, Nframes
    if  RawData==nothing
        RawData = dropmissing(CSV.read(fileName, DataFrame))
        Nreps =  Int64(RawData[:,1][end])
        Nframes = Int64(size(RawData)[1]/Nreps)  # must be same for all files
    else
        NewData = dropmissing(CSV.read(fileName, DataFrame))
        Nreps =  Nreps + Int64(NewData[:,1][end])   
        RawData = vcat(RawData, NewData)    
    end

end


# Range[i,j] is range in frame i of jth rep
Range = reshape(copy(RawData.Range), (Nframes, Nreps))
# because this is the same in every rep 
# when stalking is not stochastic (no noise in trajectory)
# just use 1 column. Also nb iRange is kluge to clip the "hovering"
# section of the approach trajectory (ie when predator gets to
# min_Δ).  Will need changing if sim parameters change.
iRange = 1:475
Range = Range[iRange,1]

# Kullback-Liebler divergences from 'true' posterior to particle distributions
KLD_particles = reshape(copy(RawData.KLD), (Nframes, Nreps))[iRange,:]  # posterior particles
KLD_sample    = reshape(copy(RawData.KLDI), (Nframes, Nreps))[iRange,:]   # sample from posterior
KLD_random    = reshape(copy(RawData.KLD0), (Nframes, Nreps))[iRange,:]     # uniform sample over mat

# posterior entropy
PosteriorEntropy = reshape(copy(RawData.PosteriorEntropy), (Nframes, Nreps))[iRange,:] 

# Posterior probability that predator is within specified range
P40 = reshape(copy(RawData.PR40), (Nframes, Nreps))[iRange,:]   # 40um
P45 = reshape(copy(RawData.PR45), (Nframes, Nreps))[iRange,:] 
P50 = reshape(copy(RawData.PR50), (Nframes, Nreps))[iRange,:]

# Proportion of posterior particles within specified range
# (posterior probability estimated by particle filter)
N40 = reshape(copy(RawData.NR40), (Nframes, Nreps))[iRange,:]   # 40um
N45 = reshape(copy(RawData.NR45), (Nframes, Nreps))[iRange,:] 
N50 = reshape(copy(RawData.NR50), (Nframes, Nreps))[iRange,:] 

# quantiles of posterior range 
# e.g. 1%  probability of predator beyond QP01 
QP01 = reshape(copy(RawData.QP01), (Nframes, Nreps))[iRange,:]  
QP05 = reshape(copy(RawData.QP05), (Nframes, Nreps))[iRange,:] 
QP25 = reshape(copy(RawData.QP25), (Nframes, Nreps))[iRange,:]  
QP50 = reshape(copy(RawData.QP50), (Nframes, Nreps))[iRange,:]  

# quantiles of particle range distribution
# e.g. 1% of particles are beyond QN01
QN01 = reshape(copy(RawData.QN01), (Nframes, Nreps))[iRange,:]  
QN05 = reshape(copy(RawData.QN05), (Nframes, Nreps))[iRange,:]  
QN25 = reshape(copy(RawData.QN25), (Nframes, Nreps))[iRange,:]  
QN50 = reshape(copy(RawData.QN50), (Nframes, Nreps))[iRange,:]  


# quantiles of posterior angle distribution
Qψ01 = reshape(copy(RawData.Qψ01), (Nframes, Nreps))[iRange,:]  
Qψ05 = reshape(copy(RawData.Qψ05), (Nframes, Nreps))[iRange,:] 
Qψ25 = reshape(copy(RawData.Qψ25), (Nframes, Nreps))[iRange,:] 
Qψ50 = reshape(copy(RawData.Qψ50), (Nframes, Nreps))[iRange,:] 
Qψ75 = reshape(copy(RawData.Qψ75), (Nframes, Nreps))[iRange,:] 
Qψ95 = reshape(copy(RawData.Qψ95), (Nframes, Nreps))[iRange,:] 
Qψ99 = reshape(copy(RawData.Qψ99), (Nframes, Nreps))[iRange,:] 

# quantiles of angles to posterior porticle locations
QΘ01 = reshape(copy(RawData.QΘ01), (Nframes, Nreps))[iRange,:]  
QΘ05 = reshape(copy(RawData.QΘ05), (Nframes, Nreps))[iRange,:] 
QΘ25 = reshape(copy(RawData.QΘ25), (Nframes, Nreps))[iRange,:] 
QΘ50 = reshape(copy(RawData.QΘ50), (Nframes, Nreps))[iRange,:] 
QΘ75 = reshape(copy(RawData.QΘ75), (Nframes, Nreps))[iRange,:] 
QΘ95 = reshape(copy(RawData.QΘ95), (Nframes, Nreps))[iRange,:] 
QΘ99 = reshape(copy(RawData.QΘ99), (Nframes, Nreps))[iRange,:] 


# M-cell (probability of target in RF)
MP = reshape(copy(RawData.MP), (Nframes, Nreps))[iRange,:]  

# M-cell particle concentration 
MN = reshape(copy(RawData.MN), (Nframes, Nreps))[iRange,:]  


#=      
# BLOCK COMMENT
# This section maps regular time sampling to regular range sampling
# using cubic splines - because data were generated per frame, but
# we are interested in inference as a function of distance to predator).
# However, in the revised 2022 version of the simulation the predator has
# constant approach speed, so this isnt necessary.

 # simulation data are saved per timestep, but we want to evaluate 
 # the observer in terms of distance to predator.  
 # convert (interpolate) to values on a regular grid of pred-prey distances
Rgrid = collect(25:.25:125)

## declare arrays to hold interpolated data
# Entropy and K-L divergence 
KLD = Array{Float64,2}(undef,length(Rgrid), Nreps)
KLD0 = Array{Float64,2}(undef,length(Rgrid), Nreps)
KLDI = Array{Float64,2}(undef,length(Rgrid), Nreps)
PENT  = Array{Float64,2}(undef,length(Rgrid), Nreps)

# probability of predator within range
PR40  = Array{Float64,2}(undef,length(Rgrid), Nreps)
PR45  = Array{Float64,2}(undef,length(Rgrid), Nreps)
PR50  = Array{Float64,2}(undef,length(Rgrid), Nreps)

NR40  = Array{Float64,2}(undef,length(Rgrid), Nreps)
NR45  = Array{Float64,2}(undef,length(Rgrid), Nreps)
NR50 = Array{Float64,2}(undef,length(Rgrid), Nreps)

# Quantiles of range
QP01  = Array{Float64,2}(undef,length(Rgrid), Nreps) # Bayesian
QP05  = Array{Float64,2}(undef,length(Rgrid), Nreps)
QP25 = Array{Float64,2}(undef,length(Rgrid), Nreps)
QP50 = Array{Float64,2}(undef,length(Rgrid), Nreps)

QN01  = Array{Float64,2}(undef,length(Rgrid), Nreps) # particle
QN05  = Array{Float64,2}(undef,length(Rgrid), Nreps)
QN25 = Array{Float64,2}(undef,length(Rgrid), Nreps)
QN50 = Array{Float64,2}(undef,length(Rgrid), Nreps)

# Quantiles of angle
Qψ01 = Array{Float64,2}(undef,length(Rgrid), Nreps)
Qψ05 = Array{Float64,2}(undef,length(Rgrid), Nreps)
Qψ25 = Array{Float64,2}(undef,length(Rgrid), Nreps)
Qψ50 = Array{Float64,2}(undef,length(Rgrid), Nreps)
Qψ75 = Array{Float64,2}(undef,length(Rgrid), Nreps)
Qψ95 = Array{Float64,2}(undef,length(Rgrid), Nreps)
Qψ99 = Array{Float64,2}(undef,length(Rgrid), Nreps)

QΘ01 = Array{Float64,2}(undef,length(Rgrid), Nreps)
QΘ05 = Array{Float64,2}(undef,length(Rgrid), Nreps)
QΘ25 = Array{Float64,2}(undef,length(Rgrid), Nreps)
QΘ50 = Array{Float64,2}(undef,length(Rgrid), Nreps)
QΘ75 = Array{Float64,2}(undef,length(Rgrid), Nreps)
QΘ95 = Array{Float64,2}(undef,length(Rgrid), Nreps)
QΘ99 = Array{Float64,2}(undef,length(Rgrid), Nreps)

# cellular proximity detector
MCP = Array{Float64,2}(undef,length(Rgrid), Nreps)
MCN = Array{Float64,2}(undef,length(Rgrid), Nreps)

for rep in 1:Nreps
    
    if minimum(R[:, rep])>25.0  # kluge for failed to go distance (rare)
        R[end,rep] = 25.0
    end

    t25 = findfirst(x -> x<=25, R[:,rep])         # time to 25um range for this rep
    i = sortperm(R[1:t25,rep])
    
    # entropy & K-L divergence
    sp = Spline1D(R[i,rep], KLD_particles[i,rep])  # interpolation function (using Dierckx) 
    KLD[:,rep] = sp(Rgrid)
    sp = Spline1D(R[i,rep], KLD_random[i,rep])  # nb -R because Dierckc requires 1st arg increasing
    KLD0[:,rep] = sp(Rgrid)
    sp = Spline1D(R[i,rep], KLD_sample[i,rep])  
    KLDI[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], Pent[i,rep])  
    PENT[:,rep] = sp(Rgrid)

    # Pr(predator in range)
    sp = Spline1D(R[i,rep], P40[i,rep])  
    PR40[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], P45[i,rep])  
    PR45[:,rep] = sp(Rgrid)
    
    sp = Spline1D(R[i,rep], P50[i,rep])  
    PR50[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], N40[i,rep])  
    NR40[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], N45[i,rep])  
    NR45[:,rep] = sp(Rgrid)
    
    sp = Spline1D(R[i,rep], N50[i,rep])  
    NR50[:,rep] = sp(Rgrid)

    # quantiles of range distn
    sp = Spline1D(R[i,rep], QP01[i,rep])  
    QP01[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], QP05[i,rep])  
    QP05[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], QP25[i,rep])  
    QP25[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], QP50[i,rep])  
    QP50[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], QN01[i,rep])  
    QN01[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], QN05[i,rep])  
    QN05[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], QN25[i,rep])  
    QN25[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], QN50[i,rep])  
    QN50[:,rep] = sp(Rgrid)

    # Quantiles of angle
    sp = Spline1D(R[i,rep], Qψ01t[i,rep])  
    Qψ01[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], Qψ05t[i,rep])  
    Qψ05[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], Qψ25t[i,rep])  
    Qψ25[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], Qψ50t[i,rep])  
    Qψ50[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], Qψ75t[i,rep])  
    Qψ75[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], Qψ95t[i,rep])  
    Qψ95[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], Qψ99t[i,rep])  
    Qψ99[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], QΘ01t[i,rep])  
    QΘ01[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], QΘ05t[i,rep])  
    QΘ05[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], QΘ25t[i,rep])  
    QΘ25[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], QΘ50t[i,rep])  
    QΘ50[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], QΘ75t[i,rep])  
    QΘ75[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], QΘ95t[i,rep])  
    QΘ95[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], QΘ99t[i,rep])  
    QΘ99[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], MCPt[i,rep])  
    MCP[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], MCNt[i,rep])  
    MCN[:,rep] = sp(Rgrid)


end    
 =#  # END BLOCK COMMENT INTERPOLATION

# indices for quantiles
iQ50 = Int64(round(Nreps/2))      # index median
iQ25 = Int64(round(Nreps/4))      # lower quartile
iQ75 = Int64(round(3*Nreps/4))    # upper quartile
iQ05 = Int64(round(0.05*Nreps))   # 5th percentile
iQ95 = Int64(round(0.95*Nreps))   # 95th percentile

# sort arrays along columns, so  iQ<q> indexes the qth quantile
# for that variable in each row (ie at each distance)
KLD_particles_sorted    = sort(KLD_particles[iRange,:], dims=2)   
KLD_random_sorted       = sort(KLD_random[iRange,:], dims=2)   
KLD_sample_sorted       = sort(KLD_sample[iRange,:], dims=2)   
PosteriorEntropy_sorted = sort(PosteriorEntropy[iRange,:], dims=2)

P40_sorted = sort(P40, dims = 2)
P45_sorted = sort(P45, dims = 2)
P50_sorted = sort(P50, dims = 2)

N40_sorted = sort(N40, dims = 2)
N45_sorted = sort(N45, dims = 2)
N50_sorted = sort(N50, dims = 2);

QP01_sorted = sort(QP01, dims = 2)
QP05_sorted = sort(QP05, dims = 2)
QP25_sorted = sort(QP25, dims = 2)
QP50_sorted = sort(QP50, dims = 2)

QN01_sorted = sort(QN01, dims = 2)
QN05_sorted = sort(QN05, dims = 2)
QN25_sorted = sort(QN25, dims = 2)
QN50_sorted = sort(QN50, dims = 2)

Qψ01_sorted = sort(Qψ01, dims = 2)
Qψ05_sorted = sort(Qψ05, dims = 2)
Qψ25_sorted = sort(Qψ25, dims = 2)
Qψ50_sorted = sort(Qψ50, dims = 2)
Qψ75_sorted = sort(Qψ75, dims = 2)
Qψ95_sorted = sort(Qψ95, dims = 2)
Qψ99_sorted = sort(Qψ99, dims = 2)

QΘ01_sorted = sort(QΘ01, dims = 2)
QΘ05_sorted = sort(QΘ05, dims = 2)
QΘ25_sorted = sort(QΘ25, dims = 2)
QΘ50_sorted = sort(QΘ50, dims = 2)
QΘ75_sorted = sort(QΘ75, dims = 2)
QΘ95_sorted = sort(QΘ95, dims = 2)
QΘ99_sorted = sort(QΘ99, dims = 2)

MP_sorted = sort(MP, dims = 2)
MN_sorted = sort(MN, dims = 2)


# meanKLD = mean(KLD, dims=2)[:,1]
# meanKLDI = mean(KLDI, dims=2)[:,1]
# meanKLD0 = mean(KL