# BayesianPlacozoan module

using Makie
using Colors
using OffsetArrays

# Param structure
# physical parameters
struct Param

  # single-cell dipole source
  ρ::Float64    # Resisitivity of seawater Ω.cm
  δ::Float64    # dipole separation cm
  I::Float64    # dipole current A

  # Johnson-Nyquist noise
  kB::Float64       # Bolzmann's constant
  T::Float64        # Kelvin
  Ω::Float64        # receptor impedance Ω
  Δf::Float64       # bandwidth
  σ::Float64        # Johnson-Nyquist noise RMS
end

# Param constructor
function Param()

  ρ = 25.0          # Resisitivity of seawater 25Ω.cm
  δ = 20.e-6*100.   # dipole separation 10μm in cm
  I = 5.0e-11*1.0e6 # dipole current 50pA. converted to μA

  # Johnson-Nyquist noise
  kB = 1.38e-23           # Bolzmann's constant
  T = 300.                # Kelvin
  Ω = 20.e6               # receptor impedance
  Δf = 1.0e3              # bandwidth
  σ = sqrt(4.0*kB*T*Δf)   # Johnson-Nyquist noise RMS

  return Param(ρ, δ, I, kB, T, Ω, Δf, σ)
end

# construct parameter structure
param = Param()

# World structure
# contains scene settings and global parameters
struct World
  nFrames::Int64
  radius::Int64
  matcolor::RGBA{Float64}
  bgcolor::RGBA{Float64}
  likelycolor::RGB{Float64}
  postcolor::RGB{Float64}
  likelysize::Array{Float64}
  postsize::Array{Float64}
  likelihood::OffsetArray     # likelihood given all receptor states
  prior::OffsetArray
  posterior::OffsetArray
  nLparticles::Int64
  nPparticles::Int64
  Lparticle::Array{Float64,2}
  Pparticle::Array{Float64,2}
  Pparticle_step::Array{Float64,2}  # particle prediction steps
  priorSD::Array{Float64}  # std. dev. of prior
  Δ::Array{Float64,1}    # closest approach of predator μm (animation parameter)
end

# World constructor
function World(nFrames::Int64, radius::Int64,
                nLparticles::Int64, nPparticles::Int64, Δ::Float64)
  indx = -radius:radius
  n_indx = length(indx)
  likelihood = OffsetArray(fill(0.0, n_indx,n_indx), indx, indx)
  prior = OffsetArray(fill(0.0, n_indx,n_indx), indx, indx)
  posterior = OffsetArray(fill(0.0, n_indx,n_indx), indx, indx)

  likelycolor = RGB(0.85, 0.65, 0.35)
  postcolor = RGB(0.5, 0.75, 1.0)
  likelysize = 4
  postsize = 4
  priorSD = 150.0

  Lparticle = fill(0.0, nLparticles,2)
  Pparticle = fill(0.0, nPparticles,2)
  Pparticle_step = zeros(nPparticles,2)
# Δ: predator is attracted to prey if it is further than this
# and repelled if it is closer; for animating stalking behaviour

  default_matcolor = RGBA(.05, .35, .35, 1.0)
  default_bgcolor  = RGBA(0.0, 0.0, 0.0, 1.0)
  return World(nFrames, radius, default_matcolor, default_bgcolor,
              likelycolor, postcolor, [likelysize], [postsize],
               likelihood, prior, posterior,
               nLparticles, nPparticles, Lparticle, Pparticle, Pparticle_step,
               [priorSD], [Δ] )
end

# electroreceptor array
# including Bayesian receptive fields (likelihoods for prey proximity)
struct Ereceptor
  N::Int64
  size::Float64  # symbol size for drawing receptors
  x::Array{Float64,1}  # receptor x-coords relative to centre of placozoan
  y::Array{Float64,1}  # receptor y-coords
  state::Array{Float64,1} # 0/1 for receptor in closed/open state
  openColor::RGB
  closedColor::RGB
  # pOpen[i] is an array containing pre-computed probability pOpen[i][j,k]
  #   that the ith receptor will be in the OPEN state if the nearest edge of
  #   the predator is at the [j,k]th grid point. This is the Bayesian
  #   receptive field of the receptor a.k.a. the likelihood function for
  #   predator proximity given that the receptor is activated (in open state).
  pOpen::Array{OffsetArray,1}
end

# Ereceptor constructor
# creates N receptors in a ring centred at (0,0)
# creates likelihood array (Bayesian receptive field) for each receptor
#    initialized to all zeros
function Ereceptor(w::World, radius::Float64, N::Int64, displaysize::Float64)

   if floor(N/4)!=N/4
     error("Number of receptors must be a multiple of 4")
   end

   # N receptors equally spaced in a ring at radius radius
   x = [radius.*(cos(2π*i/N)) for i in 1:N]
   y = [radius.*(sin(2π*i/N)) for i in 1:N]
   Open = zeros(N) # initialize receptors closed  (init. state doesn't matter)

   openColor   = RGB(1.0, 1.0, 0.25)
   closedColor = RGB(0.35, 0.45, 0.35)

   # 1d vector containing N offset arrays; ith will contain RF for ith receptor
   Lhd = Array{OffsetArray,1}(undef,N)
   indx = -w.radius:w.radius   # indices for offset array
   n_indx = length(indx)
   for i in 1:N
     Lhd[i] = OffsetArray(fill(0.0, n_indx, n_indx), indx, indx)
   end

   return Ereceptor(N, displaysize, x, y,
                    zeros(N), openColor, closedColor, Lhd)
end






# Placozoan structure
struct Placozoan
  radius::Float64
  marginwidth::Float64
  gutradius::Float64
  celldiam::Float64
  x::Array{Float64,1}  # x-ccord of centre
  y::Array{Float64,1}   # y-coord of centre
  # field[i] is pre-computed bio-electric field strength
  #   at distance i μm from edge of placozoan
  field::Array{Float64,1}
  potential::Array{Float64,1}  # in μV
  fieldrange::Int64   # number of elements in field (= max range in μm)
  receptor::Ereceptor  # electroreceptor array
  speed::Array{Float64,1}
  step::Array{Float64,1}
  color::RGBA{Float64}
  gutcolor::RGBA{Float64}
  edgecolor::RGB{Float64}
end

# placozoan constructor
# specify size and margin width
# other parameters take default values; located at origin
function Placozoan(radius::Float64, margin::Float64, Nreceptors::Int64,
                    bodycolor=RGBA(0.9, 0.75, 0.65, 0.5),
                    gutcolor = RGBA(1., 0.75, 0.75, 0.25),
                    edgecolor = RGB(0.0, 0.0, 0.0) )
    fieldrange = Int(round(radius*3))
    receptor = Ereceptor(W,radius, Nreceptors,12.0)
    return Placozoan(radius, margin, radius-margin, 12.0, [0.0], [0.0],
            fill(0.0, fieldrange), fill(0.0, fieldrange), fieldrange,
            receptor, [0.0], [0.0, 0.0],
            bodycolor, gutcolor, edgecolor )
end

# # computes Bayesian receptive fields for each receptor
# # i.e. normalized likelihood for nearest edge of predator at (x,y)
# # given that the receptor channel is open
# function precomputeBayesianRF(w::World, self::Placozoan, other::Placozoan)
#
#   for i in 1:self.receptor.N  # for each receptor
#     # precompute likelihood (open state probability) for this receptor
#     # nb likelihood of predator inside self is zero
#     # (because self must be still alive to do this computation)
#    for j in -w.radius:w.radius
#       for k in -w.radius:w.radius
#         self.receptor.pOpen[i][j,k] = sqrt(j^2+k^2) > self.radius ?
#             pOpen(sqrt((self.receptor.x[i]-j)^2 + (self.receptor.y[i]-k)^2),
#                    other.potential) : 0.0
#       end
#     end
#   end
#
# end

# computes Bayesian receptive fields for each receptor
# i.e. normalized likelihood for nearest edge of predator at (x,y)
# given that the receptor channel is open
function precomputeBayesianRF(w::World, self::Placozoan, other::Placozoan)

  # computes RFs for receptors in 1st quadrant, copies to other quadrants
  Nq = self.receptor.N ÷ 4         # receptors per quadrant
  for i in 1:Nq  # for each receptor
    # precompute likelihood (open state probability) for this receptor
    # nb likelihood of predator inside self is zero
    # (because self must be still alive to do this computation)
   for j in -w.radius:w.radius
      for k in -w.radius:w.radius

        # likelihood at (j,k)
        L = sqrt(j^2+k^2) > self.radius ?
            pOpen(sqrt((self.receptor.x[i]-j)^2 + (self.receptor.y[i]-k)^2),
                   other.potential) : 0.0
        # copy to each quadrant
        self.receptor.pOpen[i][j,k]         = L
        self.receptor.pOpen[Nq+i][-k,j]     = L
        self.receptor.pOpen[2*Nq+i][-j,-k]  = L
        self.receptor.pOpen[3*Nq+i][k,-j]   = L
      end
    end
  end

end

# function computes receptor channel Open probability
# as a function of electric field strength
# calibrated to 10% thermal noise-driven open probability for target at infinity
v0 = -param.σ*log(0.1/(1.0-0.1))
pOpenGivenFieldstrength(e) =  1.0./(1 .+ exp.(-(e.-v0)/param.σ))

# function computes single-cell dipole field strength at distance r, in μV/cm
dipoleFieldstrength(r::Float64) = 2π*param.ρ*param.I*param.δ./r.^3

# precomputes field strength and potential
# as a function of distance in μm from edge of body
# due to all dipoles in a placozoan.
# updates placozoan.field and placozoan.potential
function placozoanFieldstrength(p::Placozoan)
  for a in p.celldiam:p.celldiam:(p.gutradius - p.celldiam)
    n = round(2π*a/p.celldiam)    # number of dipoles in layer
    x = [ a*cos(2π*i/n) for i in 1:n]     # location of dipole
    y = [ a*sin(2π*i/n) for i in 1:n]
    for d in 1:p.fieldrange
      r = sqrt.(((d.+p.radius.-x).^2 + y.^2)).*1.0e-4
      p.field[d] = p.field[d] + sum(dipoleFieldstrength.(r))
    end
    # electric field in μV/cm converted to potential across 10μm receptor
    # nb 1cm = 10^4 μm
    p.potential[:] = p.field./10.0e4*10.0
  end
end

# # compute microvolts across receptor from electric field
# function microvoltsFromField(p::Placozoan)
#   V = cumsum(p.field)*1.0e-4
#   p.potential[:] = V[end].-V
# end

# electroreceptor open state probability
# as a function of distance to edge of predator
function pOpen(d, V)
   i = Int(round(d)) + 1
   if i > length(V)
     i = length(V)
   end
   return pOpenGivenFieldstrength(V[i]*1.0e-6)
 end


 # compute likelihood function in world given self.receptor state
function likelihood(w::World, self::Placozoan)

   w.likelihood .= 1.0
   for i = 1:self.receptor.N
     if self.receptor.state[i]==1
       w.likelihood .*= self.receptor.pOpen[i]
     else
       w.likelihood .*= (1.0 .- self.receptor.pOpen[i])
     end
   end

   for j in -w.radius:w.radius
     for k in -w.radius:w.radius
       if (j^2 + k^2) <= self.radius^2
         w.likelihood[j,k] = 0.0
       end
     end
   end

   w.likelihood ./= maximum(w.likelihood)

 end



 # construct sensory particles in prey margin
 # by reflectObservationg likelihood sample points through skin
 function reflectObservation(w::World, p::Placozoan)
   R = sqrt.(W.Lparticle[:,1].^2 + W.Lparticle[:,2].^2)
   r = (prey.radius .- prey.marginwidth*(R.-prey.radius)./(W.radius-prey.radius))::Array{Float64,1}
   #return (r.*xLhdSample./R, r.*yLhdSample./R)
   # observationPlot[1] = r.*W.Lparticle[:,1]./R            # update reflected sample plot
   # observationPlot[2] = r.*W.Lparticle[:,2]./R

   observation = r.*W.Lparticle./R
 end

 function reflectBelief(beliefParticle_xy)
   R = sqrt.(beliefParticle_xy[:,1].^2 + beliefParticle_xy[:,2].^2)
   r = (preyRadius .- preyMargin*(R.-preyRadius)./(matRadius-preyRadius))::Array{Float64,1}
   #return (r.*xLhdSample./R, r.*yLhdSample./R)
   beliefPlot[1] = r.*beliefParticle_xy[:,1]./R            # update reflected sample plot
   beliefPlot[2] = r.*beliefParticle_xy[:,2]./R
 end


 # Function to sample from normalized likelihood by rejection
 function sample_likelihood(w::World)

     n = 0
     while n < w.nLparticles
       candidate = rand(-w.radius:w.radius,2)
       if sqrt(candidate[1]^2 + candidate[2]^2) < w.radius
         if rand()[] < w.likelihood[candidate...]
           n = n + 1
           w.Lparticle[n, :] = candidate[:]
         end
       end
     end
 end

 function updateReceptors(prey::Placozoan, predator::Placozoan)

   # calculate receptor states
   for j = 1:length(prey.receptor.state)
      range = sqrt( (predator.x[] - prey.receptor.x[j])^2  +
                   (predator.y[] - prey.receptor.y[j])^2 ) - predator.radius

      if range < 0.0
         range = 0.0
       end

      prey.receptor.state[j] = Int(rand()[] < pOpen(range, predator.potential))
   end
 end



# rotate placozoan location through angle dψ around origin
function orbit(w::World, p::Placozoan)
  dψ = π/w.nFrames
  p.x[] =  cos(dψ)*p.x[] + sin(dψ)*p.y[]
  p.y[] = -sin(dψ)*p.x[] + cos(dψ)*p.y[]
  # C = [cos(dψ) sin(dψ); -sin(dψ) cos(dψ)]
end

# predator movement
function stalk(w::World, predator::Placozoan, prey::Placozoan)

  # predator movement
  d = sqrt(predator.x[]^2 + predator.y[]^2)  # distance from origin
  v = sign(prey.radius + predator.radius + w.Δ[] - d)#(distance between edges)-Δ.
  # pink noise motion in mat frame
  predator.step[:] = 0.8*predator.step +
        0.2*randn(2).*predator.speed[]  .+
        0.1*v*predator.speed[].*([predator.x[], predator.y[]]) ./ d

  # update predator coordinates
  predator.x[] += predator.step[1]
  predator.y[] += predator.step[2]
  orbit(w, predator)

end

# prey motion = pink noise in mat frame
#global preyStep = 0.9*preyStep + 0.1*randn(2).*preySpeed
# initialize posterior samples
# truncated Gaussian distribution of distance from mat edge
# (i.e. diffusion from edge with absorbing barrier at prey)
function initialize_posterior(w::World, p::Placozoan)

  nP = 0
  while nP < w.nPparticles
    ϕ = 2.0*π*rand(1)[]
    β = 1.0e12
    while β > (w.radius-p.radius)
      β = w.priorSD[]*abs(randn(1)[])
    end
    #candidate = -matRadius .+ sceneWidth.*rand(2)  # random point in scene
    # d = sqrt(candidate[1]^2 + candidate[2]^2) # candidate distance from origin
    # if (d>preyRadius) & (d<matRadius)
      nP = nP+1
      w.Pparticle[nP,:] =  (w.radius-β).*[cos(ϕ), sin(ϕ)]
    # end
  end
end

function posteriorPredict(w::World, predator::Placozoan)

  w.Pparticle_step .= 0.8*w.Pparticle_step +
                     1.0*randn(w.nPparticles,2).*predator.speed[]
  w.Pparticle .+= w.Pparticle_step

end

# impose boundaries on posterior particle movement
 function diffusionBoundary(w::World, prey::Placozoan)

   for j in 1:w.nPparticles

     d = sqrt(w.Pparticle[j,1]^2 + w.Pparticle[j,2]^2)

     # mat edge is a dissipative reflecting boundary
     # particles are brought to rest if they hit it
     if d>w.radius
        w.Pparticle[j,:] = w.radius.*w.Pparticle[j,:]./d
        w.Pparticle_step[j,:] = [0.0, 0.0]
     end

     # edge of prey is an absorbing barrier
     # particles are anihilated if they hit it
     # and are reborn at rest on the edge of the mat
     if d < prey.radius
        ϕ = 2.0*π*rand(1)[]
        w.Pparticle[j,:] = w.radius.*[cos(ϕ), sin(ϕ)]
        w.Pparticle_step[j,:] = [0.0, 0.0]
     end
   end
 end

 # duplicate belief particles that collide with observation particles
 function bayesBelief(w::World)

   for i in 1:w.nPparticles
     ix = Int(round(w.Pparticle[i,1]))  # x-grid coord ith particle
     for j in 1:w.nLparticles
       if ix==w.Lparticle[j,1]  # found matching x-coord
         if Int(round(w.Pparticle[i,2])) == w.Lparticle[j,2] #&y-coord
           ireplace = rand(1:w.nPparticles)[]  # pick particle to replace
           w.Pparticle[ireplace,:] = w.Pparticle[i,:] + 2.0*randn(2)
         end
       end
     end
   end
 end
