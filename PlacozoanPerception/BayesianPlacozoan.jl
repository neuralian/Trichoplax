# BayesianPlacozoan module

using Makie
using AbstractPlotting.MakieLayout
using AbstractPlotting
using Colors
using OffsetArrays
using Distributions
using ImageFiltering
using CSV
using DataFrames

const max_nLparticles = 2^14
const max_nBparticles = 2^14


# colors
# scene
colour_mat = RGBA(.05, .35, .35, 1.0)
colour_background = RGBA(68/255, 1/255, 84/255, 1.0)

# external/world particles
colour_likelihood = RGB(1.0, 0.55, 0.25)
#colour_prior = RGB(0.75, 0.45, 0.45)
#colour_posterior = RGB(0.85, 0.25, 0.25)
colour_posterior = RGB(1.00,.75,.75)

# internal/spike particles
colour_observation = :yellow



# receptors
colour_receptor_OPEN  = RGB(1.0, 1.0, 1.0)
colour_receptor_CLOSED  = RGB(0.25, 0.35, 0.25)
sizeof_receptor = 6.0

#crystal cells
vision_light = RGB(1.0, 1.0, 1.0)
vision_dark = RGB(0.0, 0.0, 0.0)
sizeof_crystal = 5.0
vision_SD = 0.8

# Particle sizes  
size_likelihood = 1.5
#size_prior = 4
size_posterior = 2.5

size_observation = 0.5
#size_prediction = 2
size_belief = 1.0

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
  I = 2.5e-12*1.0e6 # dipole current 1pA. converted to μA

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

  maxRange::Int64 # world radius = maxRange of indices for arrays
  likelihood::OffsetArray     # likelihood given all receptor states
  prior::OffsetArray
  posterior::OffsetArray
  nLparticles::Array{Int64,1}
  nBparticles::Array{Int64,1}
  Lparticle::Array{Float64,2}    # likelihood particles (samples)
  Bparticle::Array{Float64,2}   # belief (posterior) particles
  Bparticle_step::Array{Float64,2}  # particle prediction steps
  priorDensity::Array{Float64,1}
  # priormean::Float64
  # priorsd::Float64  # std. dev. of prior
  PosteriorEntropy::Array{Float64,1} # in bits, 1 per time step
  KLD::Array{Float64,1} # K-L divergence from particles to posterior
  range::Array{Float64,1} # track distance between predator & prey edges
end

# Observer constructor
function Observer(maxRange, nLparticles::Int64, nBparticles::Int64,
                 priorDensity::Float64, nFrames::Int64)

  return Observer(maxRange,
               zeros(-maxRange:maxRange, -maxRange:maxRange),
               zeros(-maxRange:maxRange, -maxRange:maxRange),
               zeros(-maxRange:maxRange, -maxRange:maxRange),
               [nLparticles], [nBparticles],
               zeros(max_nLparticles,2),
               zeros(max_nBparticles,2),
               zeros(max_nBparticles,2),
               [priorDensity],
               zeros(nFrames), zeros(nFrames), zeros(nFrames))
end

# dummy observer constructor (for constructing placozoans without observers)
function Observer()
  z1 = zeros(1)
  z2 = zeros(1,1)
  zOff = OffsetArray(z2, 0:0, 0:0)
  Observer(1, zOff, zOff, zOff, [1], [1], z2, z2, z2, [1.0], z1, z1, z1)
end




# Electroreceptor definition
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

# Electroreceptor constructor
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

# dummy electroreceptor constructor (for constructing placozoan without electroreceptors)
function Ereceptor()

  return Ereceptor(0, 0, [0], [0], zeros(1),
                 Array{OffsetArray,1}(undef,1), RGB(0,0,0), RGB(0,0,0))
end


# Crystal cell receptor definition
struct CrystalCell
  N::Int64
  size::Float64  # symbol size for drawing receptors
  x::Array{Float64,1}  # receptor x-coords relative to centre of placozoan
  y::Array{Float64,1}  # receptor y-coords
  state::Array{Float64,1} # 0/1 for receptor in closed/open state
  lineOfSight::Array{Float64,1}
  pOpenV::Array{OffsetArray,1}
  lightColor::RGB
  darkColor::RGB
end

# crystal cell constructor
function CrystalCell(worldradius::Int64, placozoancrystalmargin::Float64,
                   N::Int64, crystalSize::Float64,
                   lightColor::RGB, darkColor::RGB)
   if floor(N/4)!=N/4  ###???floor?? code for lowest?
     error("Number of receptors must be a multiple of 4")
   end

   # N receptors equally spaced in a ring at radius of (gut?/Inside margin?)
   lineOfSight = [2π*i/N for i in 0:(N-1)]  # each xtal faces radially outwards
   x = placozoancrystalmargin.*cos.(lineOfSight)
   y = placozoancrystalmargin.*sin.(lineOfSight)
   # lineOfSight = [atan(y[i],x[i]) for i in 1:N]
   # 1d vector containing N offset arrays; ith will contain RF for ith receptor
   Lhd = Array{OffsetArray,1}(undef,N)
   for i in 1:N
     Lhd[i] = zeros(-worldradius:worldradius,-worldradius:worldradius)
   end

   return CrystalCell(N, crystalSize, x, y, zeros(N), lineOfSight,
                    Lhd, vision_light, vision_dark)
end

# crystal cell dummy constructor
function CrystalCell()

  return CrystalCell(0, 0, [0], [0], zeros(1), [0],
                 Array{OffsetArray,1}(undef,1), RGB(0,0,0), RGB(0,0,0))
end

# Placozoan definition
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
  fieldrange::Int64   # number of elements in field (= max maxRange in μm)
  receptor::Ereceptor  # electroreceptor array
  photoreceptor::CrystalCell
  observer::Observer  # pointer to observer
  speed::Array{Float64,1}
  step::Array{Float64,1}
  color::RGBA{Float64}
  gutcolor::RGBA{Float64}
  edgecolor::RGB{Float64}
end

# Placozoan constructor
function Placozoan(
  radius::Int64,
  margin::Int64,
  fieldrange::Int64,
  nEreceptors::Int64,
  receptorSize::Float64,
  eRange::Int64,
  nCrystalCells::Int64,
  crystalSize::Float64,
  crystalRange::Int64,
  nLparticles,
  nBparticles,
  priorDensity::Float64,
  nFrames::Int64,
  bodycolor = RGBA(0.9, 0.75, 0.65, 0.5),
  gutcolor = RGBA(1.0, 0.65, 0.8, 0.25),
  edgecolor = RGB(0.0, 0.0, 0.0),
  )

  observer =  Observer(eRange, nLparticles, nBparticles, priorDensity, nFrames)
  
    receptor = Ereceptor( eRange, radius, nEreceptors, receptorSize,  
                          colour_receptor_OPEN, colour_receptor_CLOSED)

    crystalcell = CrystalCell(crystalRange, (radius-0.75*margin), nCrystalCells, crystalSize,
                          vision_light, vision_dark)


    if fieldrange < 1
      fieldrange = 1
    end

    return Placozoan(
      radius,
      margin,
      radius - margin,
      12.0,
      [0.0],
      [0.0],
      zeros(fieldrange),
      zeros(fieldrange),
      fieldrange,
      receptor,
      crystalcell,
      observer,
      [0.0],
      [0.0, 0.0],
      bodycolor,
      gutcolor,
      edgecolor,
    )

end # Placozoan constructor

# placozoan constructor with field but no receptors or observer
function Placozoan(radius::Int64, margin::Int64, fieldrange::Int64,
                  bodycolor::RGBA, gutcolor::RGBA, edgecolor::RGB)

   return Placozoan(radius, margin, radius-margin, 12.0, [0.0], [0.0],
     zeros(fieldrange), zeros(fieldrange), fieldrange,
     Ereceptor(), CrystalCell(), Observer(), [0.0], [0.0, 0.0],
     bodycolor, gutcolor, edgecolor )

end

function initializeObserver(p::Placozoan, nLparticles::Int64, nBparticles::Int64,
  priorDensity::Float64)

   p.observer.nLparticles[]  = nLparticles
   p.observer.nBparticles[]  = nBparticles
   p.observer.priorDensity[] = priorDensity

# +   p.observer.Lparticle = zeros(nLparticles,2)
#    p.observer.Bparticle = zeros(nBparticles,2)
#    p.observer.Bparticle_step = zeros(nBparticles,2)


   likelihood(p)           # initialize likelihood given initial receptor states
   sample_likelihood(p)    # sample from normalized likelihood
   initialize_particles(p) # draw initial sample from prior
   initialize_prior(p)     # initialize numerical Bayesian prior

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


# electroreceptor open state probability
# as a function of distance to edge of predator
function Electroreceptor_pOpen(d, V)
   i = Int(round(d)) + 1
   if i > length(V)
     i = length(V)
   end
   return pOpenGivenFieldstrength(V[i]*1.0e-6)
 end


# precompute Bayesian receptive fields for each receptor
# i.e. normalized likelihood for nearest edge of predator at (x,y)
# given that the receptor channel is open
function Ereceptor_RF(self::Placozoan, other::Placozoan)

  # computes RFs for receptors in 1st quadrant, copies to other quadrants
  Nq = self.receptor.N ÷ 4         # receptors per quadrant
  for i in 1:Nq  # for each receptor
    # precompute likelihood (open state probability) for this receptor
    # nb likelihood of predator inside self is zero
    # (because self must be still alive to do this computation)
   for j in -self.observer.maxRange:self.observer.maxRange
      for k in -self.observer.maxRange:self.observer.maxRange

        # likelihood at (j,k)
        L = sqrt(j^2+k^2) > self.radius ?
        Electroreceptor_pOpen(sqrt((self.receptor.x[i]-j)^2 + (self.receptor.y[i]-k)^2),
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


 function Vreceptor_RF(self::Placozoan)

  # computes VFs for vision/crystal cells in 1st quadrant, copies to other quadrants
  Nq = self.photoreceptor.N ÷ 4
  for i in 1:Nq  # for each receptor
   for j in -self.observer.maxRange:self.observer.maxRange
      for k in -self.observer.maxRange:self.observer.maxRange

        angleFromLineOfSight = atan(k-self.photoreceptor.y[i],j-self.photoreceptor.x[i]) -
                               2π*(i-1)/self.photoreceptor.N

        if angleFromLineOfSight > π
          angleFromLineOfSight -= 2π
        end
        if angleFromLineOfSight < -π
          angleFromLineOfSight += 2π
        end

        L = sqrt(j^2+k^2) > self.radius ? Photoreceptor_pOpen(angleFromLineOfSight) : 0.0
        # copy to each quadrant
       self.photoreceptor.pOpenV[i][j,k]         = L  #[j,k]
       self.photoreceptor.pOpenV[Nq+i][-k,j]     = L  #[-k,j]
       self.photoreceptor.pOpenV[2*Nq+i][-j,-k]  = L  #[-j,-k]
       self.photoreceptor.pOpenV[3*Nq+i][k,-j]   = L  #[k,-j]
      end
    end
  end

end


 # photoreceptor open state probability
 function Photoreceptor_pOpen(deviationFromLineofSight::Float64) #lineOfSight::Array{Float64,1},
  distribution = Normal(0, 0.8)
  peak = pdf(distribution, 0.0)
  # 50% open probability for shadow on line of sight
  p = pdf(distribution, deviationFromLineofSight)*.1/peak
  return (p)
end


 # compute likelihood given receptor states
 # option to switch off each sensory modality (Electroreception/Photoreception = false)
function likelihood(p::Placozoan, Electroreception::Bool = true, Photoreception::Bool = true)

  p.observer.likelihood .= 1.0

  if Electroreception
     for i = 1:p.receptor.N
      if p.receptor.state[i]==1
        p.observer.likelihood .*= p.receptor.pOpen[i]
      else
        p.observer.likelihood .*= (1.0 .- p.receptor.pOpen[i])
      end
    end
  end

  if Photoreception
     for i = 1:p.photoreceptor.N
      if p.photoreceptor.state[i]==1
        p.observer.likelihood .*= p.photoreceptor.pOpenV[i]
      else
        p.observer.likelihood .*= (1.0 .- p.photoreceptor.pOpenV[i])
      end
    end
  end

   for j in -p.observer.maxRange:p.observer.maxRange
     for k in -p.observer.maxRange:p.observer.maxRange
       if (j^2 + k^2) <= p.radius^2
         p.observer.likelihood[j,k] = 0.0
       end
     end
   end
  # println("m:", maximum(p.observer.likelihood))

   p.observer.likelihood ./= maximum(p.observer.likelihood)

 end


function reflect(p::Placozoan)

  # likelihood
  R = sqrt.(p.observer.Lparticle[1:p.observer.nLparticles[],1].^2 + p.observer.Lparticle[1:p.observer.nLparticles[],2].^2)
  r = (p.radius .- p.marginwidth*(R.-p.radius)./
      (p.observer.maxRange-p.radius))::Array{Float64,1}
  #return (r.*xLhdSample./R, r.*yLhdSample./R)
  # observationPlot[1] = r.*W.Lparticle[:,1]./R            # update reflected sample plot
  # observationPlot[2] = r.*W.Lparticle[:,2]./R

  observation = r.*p.observer.Lparticle[1:p.observer.nLparticles[], :]./R

  # posterior
  Rp = sqrt.(p.observer.Bparticle[1:p.observer.nBparticles[],1].^2 + p.observer.Bparticle[1:p.observer.nBparticles[],2].^2)
  rp = (p.radius .- p.marginwidth*(Rp.-p.radius)./
      (p.observer.maxRange-p.radius))::Array{Float64,1}
  #return (r.*xLhdSample./R, r.*yLhdSample./R)
  # observationPlot[1] = r.*W.Lparticle[:,1]./R            # update reflected sample plot
  # observationPlot[2] = r.*W.Lparticle[:,2]./R
  belief = rp.*p.observer.Bparticle[1:p.observer.nBparticles[],:]./Rp

  (observation, belief)

end

 # Function to sample from normalized likelihood by rejection
 function sample_likelihood(p::Placozoan)

     n = 0
     while n < p.observer.nLparticles[]
       candidate = rand(-p.observer.maxRange:p.observer.maxRange,2)
       if sqrt(candidate[1]^2 + candidate[2]^2) < p.observer.maxRange
         if rand()[] < p.observer.likelihood[candidate...]
           n = n + 1
           p.observer.Lparticle[n, :] = candidate[:]
         end
       end
     end
 end


 # electroreceptor states as function of predator location
 function electroreception(prey::Placozoan, predator::Placozoan)
   
   for j = 1:length(prey.receptor.state)
      maxRange = sqrt( (predator.x[] - prey.receptor.x[j])^2  +
                   (predator.y[] - prey.receptor.y[j])^2 ) - predator.radius

      if maxRange < 0.0
         maxRange = 0.0
       end

      prey.receptor.state[j] = Int(rand()[] < Electroreceptor_pOpen(maxRange, predator.potential))
   end
 end


 # photoreceptor states as function of predator location (approach angle)
 function photoreception(prey::Placozoan, predator::Placozoan) 

  for j = 1:length(prey.photoreceptor.state)

     angle2predator = atan( predator.y[] - prey.photoreceptor.y[j],
                               predator.x[] - prey.photoreceptor.x[j])

     deviationFromLineOfSight = angle2predator - prey.photoreceptor.lineOfSight[j]

    #println(j, ", ", prey.photoreceptor.lineOfSight[j], ", ", angle2predator, ", ", deviationFromLineOfSight)

    #  if deviationFromLineOfSight > π
    #    deviationFromLineOfSight -= 2π
    #  end
    #  if deviationFromLineOfSight < -π
    #    deviationFromLineOfSight += 2π
    #  end

     prey.photoreceptor.state[j] = Int(rand()[] < Photoreceptor_pOpen(deviationFromLineOfSight))

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

  d3 = sqrt.(prey.observer.Bparticle[1:prey.observer.nBparticles[],1].^2 + prey.observer.Bparticle[1:prey.observer.nBparticles[],2].^2)

  v3 = sign.( prey.radius  + Δ .- d3)
  prey.observer.Bparticle_step[1:prey.observer.nBparticles[],:].= 0.8*prey.observer.Bparticle_step[1:prey.observer.nBparticles[],:] +
          0.2*randn(prey.observer.nBparticles[],2).*predator.speed[]
          #  .+
          # 0.1*v3.*predator.speed[].*prey.observer.Bparticle ./ d3
  prey.observer.Bparticle[1:prey.observer.nBparticles[],:] .=  prey.observer.Bparticle[1:prey.observer.nBparticles[],:] +
                              prey.observer.Bparticle_step[1:prey.observer.nBparticles[],:]

end

function initialize_posterior_particles_Gaussian(p::Placozoan)

  nB = 0
  while nB < p.observer.nBparticles[]

    # uniform random angle + Gaussian maxRange (truncated at edge of body and mat)
    ϕ = 2.0*π*rand(1)[]
    β = 0.0
    while (β > p.observer.maxRange) | (β < p.radius)
      β = p.observer.priormean + p.observer.priorsd*randn(1)[]
    end

    nB = nB+1
    p.observer.Bparticle[nB,:] =  β.*[cos(ϕ), sin(ϕ)]
  end
end

function initialize_particles(p::Placozoan)

  p.observer.Bparticle[1:p.observer.nBparticles[],:] = samplePrior(p.observer.nBparticles[], p)


end

function bayesParticleUpdate(p::Placozoan)

  δ2 = 4.0   # squared collision maxRange
  sδ = 9.0  # scatter maxRange
  nCollision = 0
  collision = fill(0, p.observer.nBparticles[])

    # mix the posterior with the initial prior
    # this prevents the particle filter from converging fully,
    # maintains "attention" over all possible locations of predator
    # even when the posterior indicates low uncertainty about predator location.
    # (This is a known problem with particle filters - they assign zero probability 
    # density to locations where the true density is nonzero)
    nscatter = Int(round(p.observer.priorDensity[]*p.observer.nBparticles[]))
    # select particles from posterior to scatter into prior
    iscatter = rand(1:p.observer.nBparticles[], nscatter )
    p.observer.Bparticle[iscatter, :] = samplePrior(nscatter, p)
    p.observer.Bparticle_step[iscatter,:] .= 0.0



  # list Bparticles that have collided with Lparticles

     for j = 1:p.observer.nBparticles[]  # check for collisions with belief
      collided = false
      for k = 1:nCollision
        if collision[k]==j
          collided = true
        end
      end
      if !collided
        for i = 1:p.observer.nLparticles[]
        if ((p.observer.Bparticle[j, 1] - p.observer.Lparticle[i, 1])^2 +
            (p.observer.Bparticle[j, 2] - p.observer.Lparticle[i, 2])^2 < δ2)  # new collision
          nCollision = nCollision + 1
          collision[nCollision] = j
          break
        end
      end
    end
  end
  if nCollision > 0
    # each collision produces a Poisson-distributed number of new particles
    # such that expected number of new particles is p.observer.nBparticles
    newBelief = fill(0.0, p.observer.nBparticles[], 2)

    n_newparticles =
      rand(Poisson(p.observer.nBparticles[] / nCollision), nCollision)
    count = 0
     for i = 1:nCollision
       for j = 1:n_newparticles[i]
        count = count + 1
        if count <= p.observer.nBparticles[]
          R = Inf
          while R > p.observer.maxRange  # no beliefs beyond edge of world
          newBelief[count, :] =
            p.observer.Bparticle[collision[i], :] + sδ * randn(2)
            R = sqrt(newBelief[count,1]^2 + newBelief[count,2]^2)
          end
        end
      end
    end

    # kluge number of particles (normalize the discrete distribution)
    # by random particle duplication
    for i in 1:(p.observer.nBparticles[]-count)
      newBelief[count+i,:] = newBelief[rand(1:count),:]
    end

    p.observer.Bparticle[1:p.observer.nBparticles[],:] = newBelief[:,:]

  end

end


# draw n samples from prior
# uniform on mat external to prey (annulus)
function samplePrior(n, p::Placozoan)

  # uniform angles
  ϕ = 2*π*rand(n)
  # uniform across mat disc
  r =  p.radius .+  (p.observer.maxRange - p.radius).*rand(n)

  return (hcat(r.*cos.(ϕ), r.*sin.(ϕ)))

end


function initialize_prior(p::Placozoan)

  priorsum = 0.0
  p.observer.prior .= 0.0
  for i = -p.observer.maxRange:p.observer.maxRange
    for j = -p.observer.maxRange:p.observer.maxRange
      d = sqrt(i^2 + j^2)
      if ( (d>p.radius) & (d<p.observer.maxRange) )
         p.observer.prior[i, j] = 1.0
         priorsum += 1.0
      end
    end
  end
   p.observer.prior[:,:] ./= priorsum
   p.observer.posterior[:,:] = p.observer.prior[:,:]  
end

function bayesArrayUpdate(p::Placozoan)

  posteriorSum = 0.0
   for i in -p.observer.maxRange:p.observer.maxRange
     for j in -p.observer.maxRange:p.observer.maxRange
       # posterior is dynamic prior
       p.observer.posterior[i,j] *= p.observer.likelihood[i,j]
       posteriorSum += p.observer.posterior[i,j]
     end
   end
   p.observer.posterior[:,:] = 
       (1.0-p.observer.priorDensity[])*imfilter(p.observer.posterior, Kernel.gaussian(5))./posteriorSum
   p.observer.posterior[:,:] += p.observer.priorDensity.*p.observer.prior[:,:]

 end

# Utility functions

# 1-sided quantiles of range and bearing error for posterior particles 
function particleStats(p::Placozoan, predatorBearing)

  # index active particles
  N = p.observer.nBparticles[]

  # sorted squared distance from edge of prey to posterior particles
  D2 = sort(sum(p.observer.Bparticle[1:N,:].^2, dims=2), dims=1)
  # quantiles of particle distance to predator, toward prey from 1/2 (median) to 1/128 
  QD = [sqrt(D2[Int(round(N/q))]).-p.radius for q in  [2 4 20 100]]

  # bearing error for each particle
  θ = atan.(p.observer.Bparticle[1:N,2],p.observer.Bparticle[1:N,1]) .- predatorBearing

  # unwrap
  for i in 1:N
    if θ[i] > π
      θ[i] = θ[i] - 2π
    elseif    θ[i] < -π
      θ[i] = θ[i] + 2π
    end
  end

  # quantiles of particle bearing angle from predator
  θ = sort(θ)
  Qθ = hcat( [θ[N - Int(round(N/q))]*180/π for q in  [100 20 4]], [θ[Int(round(N/q))]*180/π for q in  [2 4 20 100]])

  # return quantiles + minimum distance as tuple
  #return (Q, sqrt(D2[1]).-p.radius )
  return (QD, sqrt(D2[1]).-p.radius, Qθ, θ[1]*180/π, θ[end]*180/π )

end


# plot field and receptor open state probability as a function of distance
# use CairoMakie to allow the plot to be exported to .svg or .pdf file
# save("plot.svg", plt)
function plot_sensor(p::Placozoan)

  plt = lines(p.field)
  lines!(maximum(predator.field)*
      pOpenGivenFieldstrength(predator.potential*1.0e-6))

 return(plt)
end

# Entropy of "true" posterior in bits
function entropyBits(I::Observer)
  support = findall(x->x>1.0e-6, I.posterior)
  return -sum(I.posterior[support].*log2.(I.posterior[support]))
end

#= # Entropy recorder (saves entropy on ith timestep)
function recordEntropyBits(I::Observer, i::Int64)
  I.PosteriorEntropy[i] = -sum(I.posterior.*log2.(I.posterior))
end =#

# range recorder (saves distance to predator on ith timestep)
# nb this is centre of predator from centre of prey
function recordRange(I::Observer, predator::Placozoan, i)
    I.range[i] = sqrt(predator.x[]^2 + predator.y[]^2)
end

function KLDBits_rev(I::Observer)
   KLD = 0.0
   Q = zeros(I.nBparticles[])
   sumQ = 0.0
   for k in 1:I.nBparticles[]
     i = Int64(round(I.Bparticle[k,1]))
     j = Int64(round(I.Bparticle[k,2]))
     Q[k] = I.posterior[i,j]
     sumQ += Q[k]
   end
   Q  = Q./sumQ   # normalize posterior at sample points
   for k in 1:I.nBparticles[]
     KLD = KLD - log(2,Q[k])
   end
   KLD = KLD/I.nBparticles[] - log(2,I.nBparticles[])
   return KLD
end


function KLDBits(I::Observer)
   KLD = 0.0
   Q = zeros(I.nBparticles[])
   sumQ = 0.0
   for k in 1:I.nBparticles[]
     i = max(-I.maxRange, min(I.maxRange, Int64(round(I.Bparticle[k,1]))))
     j = max(-I.maxRange, min(I.maxRange, Int64(round(I.Bparticle[k,2]))))
     Q[k] = max(I.posterior[i,j], 1.0e-6)
     sumQ += Q[k]
   end
   Q  = Q./sumQ   # normalize posterior at sample points
   KLD = sum(Q.*log2.(Q)) + log2(I.nBparticles[])
end



function recordKLDBits(I::Observer, frame::Int64)

   Q = zeros(I.nBparticles[])
   sumQ = 0.0
   for k in 1:I.nBparticles[]
    i = max(-I.maxRange, min(I.maxRange, Int64(round(I.Bparticle[k,1]))))
    j = max(-I.maxRange, min(I.maxRange, Int64(round(I.Bparticle[k,2]))))
     Q[k] = max(I.posterior[i,j], 1.0e-6)
     sumQ += Q[k]
   end
   Q  = Q./sumQ   # normalize posterior at sample points
   I.KLD[frame] = sum(Q.*log2.(Q)) + log2(I.nBparticles[])
end
