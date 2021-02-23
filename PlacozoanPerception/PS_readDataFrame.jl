
using CSV
using DataFrames
using Statistics
using Dierckx        # spline interpolation


#include("makieTheme5.jl")

# read PS data file
fileName = "PS_10_32_FEB21_ALL.csv"
D = CSV.read(fileName, DataFrame)

Nreps = Int64(D[:,1][end])

# all range vectors as columns of array
# nb Range was saved as distance between centres of predator and prey, 
#    subtract radii to get distance from prey to proximal edge of predator
R = hcat([D.Range[ (D.rep.==i),:][:,1] for i in 1:Nreps]...).-120.0.-150.0  

# K-L divergence estimates as fcns of t
KLDt = hcat([D.KLD[ (D.rep.==i),:][:,1] for i in 1:Nreps]...)
KLDIt = hcat([D.KLDI[ (D.rep.==i),:][:,1] for i in 1:Nreps]...)
KLD0t = hcat([D.KLD0[ (D.rep.==i),:][:,1] for i in 1:Nreps]...)

# posterior entropy
Pent = hcat([D.PosteriorEntropy[ (D.rep.==i),:][:,1] for i in 1:Nreps]...)

# Pr(predator within D)
PR25t = hcat([D.PR25[ (D.rep.==i),:][:,1] for i in 1:Nreps]...)
PR50t = hcat([D.PR50[ (D.rep.==i),:][:,1] for i in 1:Nreps]...)
PR100t = hcat([D.PR100[ (D.rep.==i),:][:,1] for i in 1:Nreps]...)

NR25t = hcat([D.NR25[ (D.rep.==i),:][:,1] for i in 1:Nreps]...)
NR50t = hcat([D.NR50[ (D.rep.==i),:][:,1] for i in 1:Nreps]...)
NR100t = hcat([D.NR100[ (D.rep.==i),:][:,1] for i in 1:Nreps]...)

# quantiles of predator range distribution
QP01t = hcat([D.QP01[ (D.rep.==i),:][:,1] for i in 1:Nreps]...) # Bayesian
QP05t = hcat([D.QP05[ (D.rep.==i),:][:,1] for i in 1:Nreps]...) 
QP25t = hcat([D.QP25[ (D.rep.==i),:][:,1] for i in 1:Nreps]...)
QP50t = hcat([D.QP50[ (D.rep.==i),:][:,1] for i in 1:Nreps]...)

# quantiles of predator range distribution
QN01t = hcat([D.QN01[ (D.rep.==i),:][:,1] for i in 1:Nreps]...) # particle
QN05t = hcat([D.QN05[ (D.rep.==i),:][:,1] for i in 1:Nreps]...) 
QN25t = hcat([D.QN25[ (D.rep.==i),:][:,1] for i in 1:Nreps]...)
QN50t = hcat([D.QN50[ (D.rep.==i),:][:,1] for i in 1:Nreps]...)

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
PR25  = Array{Float64,2}(undef,length(Rgrid), Nreps)
PR50  = Array{Float64,2}(undef,length(Rgrid), Nreps)
PR100 = Array{Float64,2}(undef,length(Rgrid), Nreps)

NR25  = Array{Float64,2}(undef,length(Rgrid), Nreps)
NR50  = Array{Float64,2}(undef,length(Rgrid), Nreps)
NR100 = Array{Float64,2}(undef,length(Rgrid), Nreps)

# Quantiles of range
QP01  = Array{Float64,2}(undef,length(Rgrid), Nreps) # Bayesian
QP05  = Array{Float64,2}(undef,length(Rgrid), Nreps)
QP25 = Array{Float64,2}(undef,length(Rgrid), Nreps)
QP50 = Array{Float64,2}(undef,length(Rgrid), Nreps)

QN01  = Array{Float64,2}(undef,length(Rgrid), Nreps) # particle
QN05  = Array{Float64,2}(undef,length(Rgrid), Nreps)
QN25 = Array{Float64,2}(undef,length(Rgrid), Nreps)
QN50 = Array{Float64,2}(undef,length(Rgrid), Nreps)



for rep in 1:Nreps
    
    if minimum(R[:, rep])>25.0  # kluge for failed to go distance (rare)
        R[end,rep] = 25.0
    end

    t25 = findfirst(x -> x<=25, R[:,rep])         # time to 25um range for this rep
    i = sortperm(R[1:t25,rep])
    
    # entropy & K-L divergence
    sp = Spline1D(R[i,rep], KLDt[i,rep])  # interpolation function (using Dierckx) 
    KLD[:,rep] = sp(Rgrid)
    sp = Spline1D(R[i,rep], KLD0t[i,rep])  # nb -R because Dierckc requires 1st arg increasing
    KLD0[:,rep] = sp(Rgrid)
    sp = Spline1D(R[i,rep], KLDIt[i,rep])  
    KLDI[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], Pent[i,rep])  
    PENT[:,rep] = sp(Rgrid)

    # Pr(predator in range)
    sp = Spline1D(R[i,rep], PR25t[i,rep])  
    PR25[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], PR50t[i,rep])  
    PR50[:,rep] = sp(Rgrid)
    
    sp = Spline1D(R[i,rep], PR100t[i,rep])  
    PR100[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], NR25t[i,rep])  
    NR25[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], NR50t[i,rep])  
    NR50[:,rep] = sp(Rgrid)
    
    sp = Spline1D(R[i,rep], NR100t[i,rep])  
    NR100[:,rep] = sp(Rgrid)

    # quantiles of range distn
    sp = Spline1D(R[i,rep], QP01t[i,rep])  
    QP01[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], QP05t[i,rep])  
    QP05[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], QP25t[i,rep])  
    QP25[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], QP50t[i,rep])  
    QP50[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], QN01t[i,rep])  
    QN01[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], QN05t[i,rep])  
    QN05[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], QN25t[i,rep])  
    QN25[:,rep] = sp(Rgrid)

    sp = Spline1D(R[i,rep], QN50t[i,rep])  
    QN50[:,rep] = sp(Rgrid)


end    


# plot quantiles
iQ50 = Int64(round(Nreps/2))      # index median
iQ25 = Int64(round(Nreps/4))      # lower quartile
iQ75 = Int64(round(3*Nreps/4))    # upper quartile
iQ05  = Int64(round(0.05*Nreps))   # 5th percentile
iQ95 = Int64(round(0.95*Nreps))   # 95th percentile


KLDS = sort(KLD, dims=2)   # sort rows
KLD0S = sort(KLD0, dims=2)   # sort rows
KLDIS = sort(KLDI, dims=2)   # sort rows
PENTS = sort(PENT, dims=2)

PR25S = sort(PR25, dims = 2)
PR50S = sort(PR50, dims = 2)
PR100S = sort(PR100, dims = 2)

NR25S = sort(NR25, dims = 2)
NR50S = sort(NR50, dims = 2)
NR100S = sort(NR100, dims = 2);

QP01S = sort(QP01, dims = 2)
QP05S = sort(QP05, dims = 2)
QP25S = sort(QP25, dims = 2)
QP50S = sort(QP50, dims = 2)

QN01S = sort(QN01, dims = 2)
QN05S = sort(QN05, dims = 2)
QN25S = sort(QN25, dims = 2)
QN50S = sort(QN50, dims = 2)


# meanKLD = mean(KLD, dims=2)[:,1]
# meanKLDI = mean(KLDI, dims=2)[:,1]
# meanKLD0 = mean(KL