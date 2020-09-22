# BayesianPlacozoan module

using Makie
using Colors
using OffsetArrays
using Distributions

# colors
# scene
colour_mat = RGBA(.05, .35, .35, 1.0)
colour_background = RGBA(0.0, 0.0, 0.0, 1.0)

# external/world particles
colour_likelihood = RGB(1.0, 0.55, 0.25)
#colour_prior = RGB(0.75, 0.45, 0.45)
#colour_posterior = RGB(0.85, 0.25, 0.25)
colour_posterior = RGB(.75,.75,1.0)

# internal/spike particles
colour_observation = :yellow


# receptors
colour_receptor_OPEN  = RGB(1.0, 1.0, 0.25)
colour_receptor_CLOSED  = RGB(0.35, 0.45, 0.35)
sizeof_receptor = 8.0

# Particle sizes
size_likelihood = 2.0
#size_prior = 4
size_posterior = 2.5

size_observation = 1.5
#size_prediction = 2
size_belief = 2.0

# Physics structure
# contains physical parameters
struct Physics

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

# Physics constructor
function Physics()

  ρ = 25.0          # Resisitivity of seawater 25Ω.cm
  δ = 20.e-6*100.   # dipole separation 10μm in cm
  I = 1.0e-12*1.0e6 # dipole current 1pA. converted to μA

  # Johnson-Nyquist noise
  kB = 1.38e-23           # Bolzmann's constant
  T = 300.                # Kelvin
  Ω = 20.e6               # receptor impedance
  Δf = 1.0e3              # bandwidth
  σ = sqrt(4.0*kB*T*Δf)   # Johnson-Nyquist noise RMS

  return Physics(ρ, δ, I, kB, T, Ω, Δf, σ)

end

# construct physics
physics = Physics()


struct Observer

  range::Int64 # world radius = range of indices for arrays
  likelihood::OffsetArray     # likelihood given all receptor states
  posterior::OffsetArray
  nLparticles::Int64
  nBparticles::Int64
  Lparticle::Array{Float64,2}    # likelihood particles (samples)
  Bparticle::Array{Float64,2}   # belief (posterior) particles
  Bparticle_step::Array{Float64,2}  # particle prediction steps
  priorSD::Float64  # std. dev. of prior

end

# Observer constructor
function Observer(range,
                  nLparticles::Int64, nBparticles::Int64, priorSD::Float64)

  likelihood = zeros(-range:range, -range:range)
  posterior = zeros(-range:range, -range:range)

  Lparticle = zeros(nLparticles,2)
  Bparticle = zeros(nBparticles,2)
  Bparticle_step = zeros(nBparticles,2)

  return Observer(range, likelihood, posterior, nLparticles, nBparticles,
               Lparticle, Bparticle, Bparticle_step, priorSD)
end

# dummy observer constructor
# (for constructing placozoans without observers)
function Observer()
  z = zeros(1,1)
  zOff = OffsetArray(z, 0:0, 0:0)
  Observer(1, zOff, zOff, 1, 1, z, z, z, 1.0)
end


# electroreceptor array
# including Bayesian receptive fields (likelihoods for prey proximity)
struct Ereceptor
  N::Int64
  size::Float64  # symbol size for drawing receptors
  x::Array{Float64,1}  # receptor x-coords relative to centre of placozoan
  y::Array{Float64,1}  # receptor y-coords
  state::Array{Float64,1} # 0/1 for receptor in closed/open state

  # pOpen[i] is an array containing pre-computed probability pOpen[i][j,k]
  #   that the ith receptor will be in the OPEN state if the nearest edge of
  #   the predator is at the [j,k]th grid point. This is the Bayesian
  #   receptive field of the receptor a.k.a. the likelihood function for
  #   predator proximity given that the receptor is ON
  # (and no other observations are available)
  pOpen::Array{OffsetArray,1}
  openColor::RGB
  closedColor::RGB
end

# Ereceptor constructor
# creates N receptors in a ring centred at (0,0)
# creates likelihood array (Bayesian receptive field) for each receptor
#    initialized to all zeros
function Ereceptor(worldradius::Int64, placozoanradius::Int64,
                   N::Int64, receptorSize::Float64,
                   openColor::RGB, closedColor::RGB)

   if floor(N/4)!=N/4
     error("Number of receptors must be a multiple of 4")
   end

   # N receptors equally spaced in a ring at radius radius
   x = [placozoanradius.*(cos(2π*i/N)) for i in 1:N]
   y = [placozoanradius.*(sin(2π*i/N)) for i in 1:N]

   # 1d vector containing N offset arrays; ith will contain RF for ith receptor
   Lhd = Array{OffsetArray,1}(undef,N)
   for i in 1:N
     Lhd[i] = zeros(-worldradius:worldradius,-worldradius:worldradius)
   end

   return Ereceptor(N, receptorSize, x, y, zeros(N),
                    Lhd, colour_receptor_OPEN, colour_receptor_CLOSED)
end

# dummy Ereceptor constructor
# for constructing placozoan without receptors
function Ereceptor()

  return Ereceptor(0, 0, [0], [0], zeros(1),
                 Array{OffsetArray,1}(undef,1), RGB(0,0,0), RGB(0,0,0))
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
  observer::Observer
  speed::Array{Float64,1}
  step::Array{Float64,1}
  color::RGBA{Float64}
  gutcolor::RGBA{Float64}
  edgecolor::RGB{Float64}
end

# placozoan constructor
#   radius::Int64 = placozoan disc radius
#   margin::Int64 = margin width (observer zone)
#   nEreceptors::Int64 = number of electroreceptors
#   eRange::Int64 = range of observations (typically the radius of the 'world')
#
function Placozoan(radius::Int64, margin::Int64, fieldrange::Int64,
                    nEreceptors::Int64, receptorSize::Float64, eRange::Int64,
                    nLparticles, nBparticles, priorSD::Float64,
                    bodycolor=RGBA(0.9, 0.75, 0.65, 0.5),
                    gutcolor = RGBA(1., 0.65, 0.8, 0.25),
                    edgecolor = RGB(0.0, 0.0, 0.0) )

    observer = Observer(eRange, nLparticles, nBparticles, priorSD)
    receptor = Ereceptor(eRange,radius, Nreceptors, receptorSize,
              colour_receptor_OPEN, colour_receptor_CLOSED)

    if fieldrange<1 fieldrange = 1; end

    return Placozoan(radius, margin, radius-margin, 12.0, [0.0], [0.0],
            zeros(fieldrange), zeros(fieldrange), fieldrange,
            receptor, observer, [0.0], [0.0, 0.0],
            bodycolor, gutcolor, edgecolor )
end

# placozoan constructor with field but no receptors or observer
function Placozoan(radius::Int64, margin::Int64, fieldrange::Int64,
                  bodycolor::RGBA, gutcolor::RGBA, edgecolor::RGB)

   return Placozoan(radius, margin, radius-margin, 12.0, [0.0], [0.0],
     zeros(fieldrange), zeros(fieldrange), fieldrange,
     Ereceptor(), Observer(), [0.0], [0.0, 0.0],
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
function precomputeBayesianRF(self::Placozoan, other::Placozoan)

  # computes RFs for receptors in 1st quadrant, copies to other quadrants
  Nq = self.receptor.N ÷ 4         # receptors per quadrant
  for i in 1:Nq  # for each receptor
    # precompute likelihood (open state probability) for this receptor
    # nb likelihood of predator inside self is zero
    # (because self must be still alive to do this computation)
   for j in -self.observer.range:self.observer.range
      for k in -self.observer.range:self.observer.range

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
v0 = -physics.σ*log(0.1/(1.0-0.1))
pOpenGivenFieldstrength(e) =  1.0./(1 .+ exp.(-(e.-v0)/physics.σ))

# function computes single-cell dipole field strength at distance r, in μV/cm
dipoleFieldstrength(r::Float64) = 2π*physics.ρ*physics.I*physics.δ./r.^3

# precomputes field strength and potential
# as a function of distance in μm from edge of body
# due to all dipoles in a placozoan.
# updates placozoan.field and placozoan.potential
function placozoanFieldstrength!(p::Placozoan)
  for a in p.celldiam:p.celldiam:(p.gutradius - p.celldiam)
    n = round(2π*a/p.celldiam)    # number of dipoles in layer
    x = [ a*cos(2π*i/n) for i in 1:n]     # location of dipole
    y = [ a*sin(2π*i/n) for i in 1:n]
    for d in 1:p.fieldrange
      r = sqrt.(((d.+p.radius.-x).^2 + y.^2)).*1.0e-4
      #r = sqrt.(((d.-x).^2 + y.^2)).*1.0e-4
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


 # compute likelihood given receptor states
function likelihood(p::Placozoan)

   p.observer.likelihood .= 1.0
   for i = 1:p.receptor.N
     if p.receptor.state[i]==1
       p.observer.likelihood .*= p.receptor.pOpen[i]
     else
       p.observer.likelihood .*= (1.0 .- p.receptor.pOpen[i])
     end
   end

   for j in -p.observer.range:p.observer.range
     for k in -p.observer.range:p.observer.range
       if (j^2 + k^2) <= p.radius^2
         p.observer.likelihood[j,k] = 0.0
       end
     end
   end

   p.observer.likelihood ./= maximum(p.observer.likelihood)

 end


function reflect(p::Placozoan)

  # likelihood
  R = sqrt.(p.observer.Lparticle[:,1].^2 + p.observer.Lparticle[:,2].^2)
  r = (p.radius .- p.marginwidth*(R.-p.radius)./
      (p.observer.range-prey.radius))::Array{Float64,1}
  #return (r.*xLhdSample./R, r.*yLhdSample./R)
  # observationPlot[1] = r.*W.Lparticle[:,1]./R            # update reflected sample plot
  # observationPlot[2] = r.*W.Lparticle[:,2]./R

  observation = r.*p.observer.Lparticle./R

  # posterior
  Rp = sqrt.(p.observer.Bparticle[:,1].^2 + p.observer.Bparticle[:,2].^2)
  rp = (p.radius .- p.marginwidth*(Rp.-p.radius)./
      (p.observer.range-p.radius))::Array{Float64,1}
  #return (r.*xLhdSample./R, r.*yLhdSample./R)
  # observationPlot[1] = r.*W.Lparticle[:,1]./R            # update reflected sample plot
  # observationPlot[2] = r.*W.Lparticle[:,2]./R
  belief = rp.*p.observer.Bparticle./Rp

  (observation, belief)

end

 # Function to sample from normalized likelihood by rejection
 function sample_likelihood(p::Placozoan)

     n = 0
     while n < p.observer.nLparticles
       candidate = rand(-p.observer.range:p.observer.range,2)
       if sqrt(candidate[1]^2 + candidate[2]^2) < p.observer.range
         if rand()[] < p.observer.likelihood[candidate...]
           n = n + 1
           p.observer.Lparticle[n, :] = candidate[:]
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
function orbit(dψ::Float64, p::Placozoan)
  p.x[] =  cos(dψ)*p.x[] + sin(dψ)*p.y[]
  p.y[] = -sin(dψ)*p.x[] + cos(dψ)*p.y[]
  # C = [cos(dψ) sin(dψ); -sin(dψ) cos(dψ)]
end

# rotate particles through angle dψ around origin
function orbit(dψ::Float64, p::Array{Float64,2})
  # p.x[] =  cos(dψ)*p.x[] + sin(dψ)*p.y[]
  # p.y[] = -sin(dψ)*p.x[] + cos(dψ)*p.y[]
  p = p*[cos(dψ) -sin(dψ); sin(dψ) cos(dψ)]
   nothing
end

# predator movement
function stalk(predator::Placozoan, prey::Placozoan, Δ::Float64)

  # predator movement
  d = sqrt(predator.x[]^2 + predator.y[]^2)  # distance from origin
  v = sign(prey.radius + predator.radius + Δ - d)#(distance between edges)-Δ.
  # pink noise motion in mat frame
  predator.step[:] = 0.8*predator.step +
                    0.2*randn(2).*predator.speed[]  .+
                    0.1*v*predator.speed[].*([predator.x[], predator.y[]]) ./ d

  # update predator coordinates
  predator.x[] += predator.step[1]
  predator.y[] += predator.step[2]
  #orbit(π/w.nFrames, predator)

  # d2 = sqrt.(prey.observer.Pparticle[:,1].^2 + prey.observer.Pparticle[:,2].^2)
  # v2 = sign.( prey.radius  + Δ .- d2)
  # prey.observer.Pparticle_step .= 0.8*prey.observer.Pparticle_step +
  #               0.2*randn(prey.observer.nPparticles,2).*predator.speed[] .+
  #             0.1*v2.*predator.speed[].*prey.observer.Pparticle ./ d2
  # prey.observer.Pparticle .=  prey.observer.Pparticle +
  #                             prey.observer.Pparticle_step
  #orbit(π/w.nFrames, w.Pparticle)

  d3 = sqrt.(prey.observer.Bparticle[:,1].^2 + prey.observer.Bparticle[:,2].^2)
  # for i in 1:w.nBparticles
  #    if d3[i] > w.radius
  #     w.Bparticle[i,:] ./ d3[i]
  #     d3[i] = w.radius
  #     w.Bparticle_step[i,:] = [0.0, 0.0]
  #   end
  # end
  v3 = sign.( prey.radius  + Δ .- d3)
  prey.observer.Bparticle_step .= 0.8*prey.observer.Bparticle_step +
          0.2*randn(prey.observer.nBparticles[],2).*predator.speed[] .+
          0.1*v3.*predator.speed[].*prey.observer.Bparticle ./ d3
  prey.observer.Bparticle .=  prey.observer.Bparticle +
                              prey.observer.Bparticle_step

end

# # initialize posterior samples
# # truncated Gaussian distribution of distance from mat edge
# # (i.e. diffusion from edge with absorbing barrier at prey)
# function initialize_prior_Gaussian(p::Placozoan)
#
#   nP = 0
#   while nP < p.observer.nPparticles
#     ϕ = 2.0*π*rand(1)[]
#     β = 1.0e12
#     while β > (p.observer.range - p.radius)
#       β = p.observer.priorSD[]*abs(randn(1)[])
#     end
#     #candidate = -matRadius .+ sceneWidth.*rand(2)  # random point in scene
#     # d = sqrt(candidate[1]^2 + candidate[2]^2) # candidate distance from origin
#     # if (d>preyRadius) & (d<matRadius)
#       nP = nP+1
#       p.observer.Pparticle[nP,:] =  (p.observer.range-β).*[cos(ϕ), sin(ϕ)]
#     # end
#   end
# end

function initialize_posterior_Gaussian(p::Placozoan)

  nB = 0
  while nB < p.observer.nBparticles
    ϕ = 2.0*π*rand(1)[]
    β = 1.0e12
    while β > (p.observer.range - p.radius)
      β = p.observer.priorSD[]*abs(randn(1)[])
    end
    #candidate = -matRadius .+ sceneWidth.*rand(2)  # random point in scene
    # d = sqrt(candidate[1]^2 + candidate[2]^2) # candidate distance from origin
    # if (d>preyRadius) & (d<matRadius)
      nB = nB+1
      p.observer.Bparticle[nB,:] =  (p.observer.range - β).*[cos(ϕ), sin(ϕ)]
    # end
  end
end

# function initialize_prior_uniform(p::Placozoan)
#
#     ϕ = 2.0*π*rand(p.observer.nPparticles)
#     β = p.radius .+ rand(p.observer.nPparticles).*(p.observer.range - p.radius)
#     p.observer.Pparticle[:] =  hcat(β.*cos.(ϕ), β.*sin.(ϕ))
#
# end

# function posteriorPredict(w::World, predator::Placozoan)
#
#   w.Pparticle_step .= 0.8*w.Pparticle_step +
#                      1.0*randn(w.nPparticles,2).*predator.speed[]
#   w.Pparticle .=  w.Pparticle + w.Pparticle_step
#   #orbit(π/w.nFrames, w.Pparticle)
#
# end

# # impose boundaries on posterior particle movement
#  function steadyPrior(p::Placozoan)
#
#      pDie = 0.025
#    for j in 1:p.observer.nPparticles
#
#      d = sqrt(p.observer.Pparticle[j,1]^2 + p.observer.Pparticle[j,2]^2)
#
#      # mat edge is a dissipative reflecting boundary
#      # particles are brought to rest if they hit it
#      if d>p.observer.range
#         p.observer.Pparticle[j,:] =
#                          p.observer.range.*p.observer.Pparticle[j,:]./d
#         p.observer.Pparticle_step[j,:] = [0.0, 0.0]
#
#      # edge of prey is an absorbing barrier
#      # particles are anihilated if they hit it
#      # and are reborn at rest on the edge of the mat
#    elseif d < prey.radius
#         ϕ = 2.0*π*rand(1)[]
#         p.observer.Pparticle[j,:] = p.observer.range.*[cos(ϕ), sin(ϕ)]
#         p.observer.Pparticle_step[j,:] = [0.0, 0.0]
#
#      # particles die at random
#      # and are replaced by new particles at the mat edge
#    elseif rand()[]<pDie
#          ϕ = 2*π*rand()[]
#          p.observer.Pparticle[j,:] =
#               [p.observer.range*cos(ϕ), p.observer.range*sin(ϕ)]
#          p.observer.Pparticle_step[j,:] = [0.0, 0.0]
#       end
#
#    end
#  end

# Bayes update rule: Create a new belief particle
# when an observation particle collides with a hypothesis particle
# A hypothesis is a belief or a prior belief
# function bayesUpdate(p::Placozoan)
#
#   δ2 = 1.0
#   for i in 1:p.observer.nLparticles
#
#     for j in 1:p.observer.nPparticles  # check for collisions with prior
#       if (p.observer.Pparticle[j,1] - p.observer.Lparticle[i,1])^2 +
#          (p.observer.Pparticle[j,2] - p.observer.Lparticle[i,2])^2 < δ2   # collision
#             ireplace = rand(1:p.observer.nBparticles[])[]  # pick particle to replace
#             p.observer.Bparticle[ireplace,:] =
#                                  p.observer.Lparticle[i,:] + 8.0*randn(2)
#       end
#     end
#
#       for j in 1:p.observer.nBparticles[]  # check for collisions with belief
#         if (p.observer.Bparticle[j,1] - p.observer.Lparticle[i,1])^2 +
#            (p.observer.Bparticle[j,2] - p.observer.Lparticle[i,2])^2 < δ2   # collision
#               ireplace = rand(1:p.observer.nBparticles[])[]  # pick particle to replace
#               p.observer.Bparticle[ireplace,:] =
#                    p.observer.Lparticle[i,:] + 8.0*randn(2)
#
#         else # if no collision then kill
#           if rand()[] < 1.0e-6
#           ϕ = 2*π*rand()[]
#           p.observer.Bparticle[j,:] =
#                [p.observer.range*cos(ϕ), p.observer.range*sin(ϕ)]
#           p.observer.Bparticle_step[j,:] = [0.0, 0.0]
#        end
#         end
#
#       end   # for Bparticles
#
#     end  # for Lparticles
#
#   end

  function bayesUpdate(p::Placozoan)

    δ2 = 9.0   # squared collision range
    nCollision = 0
    collision = fill(0, 10*p.observer.nLparticles)
    # list Bparticles that have collided with Lparticles
    for i in 1:p.observer.nLparticles
      for j in 1:p.observer.nBparticles[]  # check for collisions with belief
        if (p.observer.Bparticle[j,1] - p.observer.Lparticle[i,1])^2 +
          (p.observer.Bparticle[j,2] - p.observer.Lparticle[i,2])^2 < δ2   # collision
           nCollision = nCollision + 1
           collision[nCollision] = j
         end
       end
     end
     if nCollision > 0
     # each collision produces a Poisson-distributed number of new particles
     # such that expected number of new particles is p.observer.nBparticles
     newBelief = fill(0.0, p.observer.nBparticles, 2)

     n_newparticles = rand(Poisson(p.observer.nBparticles/nCollision), nCollision)
     count = 0
     for i in 1:nCollision
        for j in 1:n_newparticles[i]
          count = count + 1
          if count <= p.observer.nBparticles
          newBelief[count,:] = p.observer.Bparticle[collision[i],:] + 12.0*randn(2)
          end
        end
      end

    p.observer.Bparticle[:] = newBelief[:]

    # scatter 100S% of Bparticles on boundary
    S = 0.001
    nscatter = Int(round(S*p.observer.nBparticles))
    iscatter = rand(1:p.observer.nBparticles, nscatter )
    ϕ = 2*π*rand(nscatter)
    r =  p.observer.range .-  rand(Exponential(50.0),nscatter)
    for i in 1:nscatter
      p.observer.Bparticle[iscatter[i], :] = [r[i]*cos(ϕ[i]), r[i]*sin(ϕ[i])]
      p.observer.Bparticle_step[iscatter[i],:] = [0.0, 0.0]
    end

end

end
