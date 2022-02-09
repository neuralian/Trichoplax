# BayesianPlacozoan module
# types and functions &c used by PlacozoanStalker.jl

using GLMakie
#using AbstractPlotting.MakieLayout
#using AbstractPlotting
using Colors
using OffsetArrays
using Distributions
using Random
using ImageFiltering
using CSV
using DataFrames
using PolygonOps    

diffuseCoef = 3.0


const max_nLparticles = 2^14
const max_nPparticles = 2^14


# colors
# scene
const colour_mat =  RGB(.1, .45, .4) # "#87966c"
const colour_background = RGB(0.1, 0.1, 0.1)
const title_color = RGB(.5, .5, .75)

# external/world particles
const colour_likelihood = "#f5cd4c" #RGB(1.0,.85, 0.65)
#colour_prior = RGB(0.75, 0.45, 0.45)
#colour_posterior = RGB(0.85, 0.25, 0.25)
const colour_posterior = RGB(.75,.25, 0.25) #"#a02c7d" #

# internal/spike particles
const colour_observation = :yellow

# placozoans
const gutcolor = RGB(0.25, 0.25, 0.25)



# receptors
const colour_receptor_OPEN  = RGB(255/255,235/255,50/255)
const colour_receptor_CLOSED  = RGB(100/255,120/255,75/255)
const sizeof_receptor = 6.0

#crystal cells
const vision_light = RGB(1.0, 1.0, 1.0)
const vision_dark = RGB(0.0, 0.0, 0.0)
const sizeof_crystal = 5.0
const vision_SD = 0.8

# "Mauthner" cell
const mcell_radius = 2.5
const mcell_inset = 2.0*mcell_radius  # inset of M-cell centre from edge of animal

# Particle sizes  
const size_likelihood = 2
#size_prior = 4
const size_posterior = 3

const size_observation = 0.5
#size_prediction = 2
const size_belief = 2.0

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


# "Mauthner" cell soma
# estimates probability of imminent threat
struct Mcell
  d::Float64           # distance from Placozoan centre to M-cell centre 
  r::Float64           # cell radius
  x::Array{Float64,1}         # x-y coords change
  y::Array{Float64,1}
end


struct Observer

  minRange::Int64  # placozoan radius (distance/index to proximal edge of mat)
  maxRange::Int64  # mat radius (distance/index to edge of mat/world)
  N::Int64         # number of sample points on mat
  log2N::Float64   # log2 of N
  likelihood::OffsetArray     # likelihood given all receptor states
  prior::OffsetArray
  posterior::OffsetArray
  nLparticles::Array{Int64,1}
  nPparticles::Array{Int64,1}
  Lparticle::Array{Float64,2}    # likelihood particles 
  Pparticle::Array{Float64,2}    # posterior particles
  Sparticle::Array{Float64,2}    # sensory  particles (reflected likelihood particles)
  Bparticle::Array{Float64,2}    # belief particles (reflected posterior particles)
  Pparticle_step::Array{Float64,2}  # particle prediction steps
  posteriorDeathRate::Array{Float64,1}   # probability that a posterior particle will die and be replaced by sample from prior per frame
  diffuseCoef::Array{Float64,1}          # posterior diffusion coefficient
  collisionRadius::Array{Float64,1}
  burnIn::Int64
  # priormean::Float64
  # priorsd::Float64  # std. dev. of prior
  PosteriorEntropy::Array{Float64,1} # in bits, 1 per time step
  KLD::Array{Float64,1}   # K-L divergence from particles to posterior
  KLD0::Array{Float64,1}  # K-L divergence from uniform sample to posterior
  KLDI::Array{Float64,1}  # K-L divergence of random sample from posterior
  Δ::Array{Float64,1} # track distance between predator & prey edges
end

# Observer constructor
function Observer(minRange, maxRange, nLparticles::Int64, nPparticles::Int64,
                  posteriorDeathRate::Float64, diffuseCoef::Float64, collisionRadius::Float64,  nFrames::Int64)


  # count sample points on mat
  N = 0
  for i in -maxRange:maxRange
    for j in -maxRange:maxRange
      d = sqrt(i^2 + j^2)
      if (d >= minRange) & (d <= maxRange)
        N = N+1
      end
    end
  end



  return Observer(minRange, maxRange, N, log2(N),
               zeros(-maxRange:maxRange, -maxRange:maxRange),
               zeros(-maxRange:maxRange, -maxRange:maxRange),
               zeros(-maxRange:maxRange, -maxRange:maxRange),
               [nLparticles], [nPparticles],
               zeros(nPparticles,2),
               zeros(nPparticles,2),
               zeros(nPparticles,2),
               zeros(nPparticles,2),
               zeros(nPparticles,2),
               [posteriorDeathRate],
               [diffuseCoef],
               [collisionRadius],
               32,
               zeros(nFrames), zeros(nFrames), zeros(nFrames),
               zeros(nFrames), zeros(nFrames))
end

# dummy observer constructor (for constructing placozoans without observers)
function Observer()
  z1 = zeros(1)
  z2 = zeros(1,1)
  zOff = OffsetArray(z2, 0:0, 0:0)
  Observer(1, 1, 1, 1.0, zOff, zOff, zOff, [1], [1], z2, z2, z2, z2, z2, 
          [0.0],[0.0], [0.0], 0, z1, z1, z1, z1, z1)
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
  observer::Observer  # Bayesian particle filter
  mcell::Mcell        # "Mauthner" neuron
  speed::Array{Float64,1}
  ω::Array{Float64,1}   # omega, angular velocity
  step::Array{Float64,1}
  color::RGBA{Float64}
  gutcolor::RGBA{Float64}
  edgecolor::RGB{Float64}
end

# Placozoan constructor
# With observer
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
  nPparticles,
  posteriorDeathRate::Float64,
  diffuseCoef::Float64,
  collisionRadius::Float64,
  nFrames::Int64,
  bodycolor = RGBA(0.9, 0.75, 0.65, 0.5),
  gutcolor = gutcolor,
  edgecolor = RGB(0.0, 0.0, 0.0),
  )

  observer =  Observer(radius, eRange, nLparticles, nPparticles, 
                        posteriorDeathRate, diffuseCoef, collisionRadius, nFrames)
  
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
      Mcell(radius-mcell_inset, mcell_radius, [0.0], [0.0]),
      [0.0],  
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
     Ereceptor(), CrystalCell(), Observer(), Mcell(0.0, 0.0, [0.0], [0.0]), [0.0], [0.0], [0.0, 0.0],
     bodycolor, gutcolor, edgecolor )

end

function initializeObserver(p::Placozoan, nLparticles::Int64, nPparticles::Int64,
  priorDensity::Float64)

   p.observer.nLparticles[]  = nLparticles
   p.observer.nPparticles[]  = nPparticles
   p.observer.priorDensity[] = priorDensity

# +   p.observer.Lparticle = zeros(nLparticles,2)
#    p.observer.Pparticle = zeros(nPparticles,2)
#    p.observer.Pparticle_step = zeros(nPparticles,2)


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
function likelihood(p::Placozoan, Electroreception::Bool = true, Photoreception::Bool = false)

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
       r = sqrt(j^2 + k^2)
       if  (r<= p.radius) | (r>=p.observer.maxRange)
         p.observer.likelihood[j,k] = 0.0
       end
     end
   end
  # println("m:", maximum(p.observer.likelihood))

   p.observer.likelihood ./= maximum(p.observer.likelihood)

 end


# map likelihood particles and posterior particles from mat onto the marginal zone of the placozoan
function reflectParticles!(p::Placozoan)

  # likelihood particles
  R = sqrt.(p.observer.Lparticle[1:p.observer.nLparticles[],1].^2 + p.observer.Lparticle[1:p.observer.nLparticles[],2].^2)
  r = (p.radius .- p.marginwidth*(R.-p.radius)./
      (p.observer.maxRange-p.radius))::Array{Float64,1}
  p.observer.Sparticle[1:p.observer.nLparticles[],:] = r.*p.observer.Lparticle[1:p.observer.nLparticles[], :]./R



  # posterior particles
  Rp = sqrt.(p.observer.Pparticle[1:p.observer.nPparticles[],1].^2 + p.observer.Pparticle[1:p.observer.nPparticles[],2].^2)
  rp = (p.radius .- p.marginwidth*(Rp.-p.radius)./ (p.observer.maxRange-p.radius))::Array{Float64,1}
  p.observer.Bparticle[1:p.observer.nPparticles[],:] = rp.*p.observer.Pparticle[1:p.observer.nPparticles[],:]./Rp


end

# map likelihood and posterior density from mat onto the marginal zone of the placozoan
function reflectArrays!(p::Placozoan)

  # Likelihood
  @inbounds for i in -p.radius:p.radius
    @inbounds for j in -p.radius:p.radius
      r = sqrt(i^2+j^2)
      if (r<p.radius) & (r > (p.radius-p.marginwidth))  # in marginal zone
        # project map location to real-world location
       # R = p.radius + (p.radius - r*(p.observer.maxRange - p.radius))/p.marginwidth 

        R = (p.radius - r)*(p.observer.maxRange-p.radius)/ p.marginwidth + p.radius
                
        #println(i, ", ", j, ", ", r, ", ", R, ", ", R/r)
        iproj = Int64(round(i*R/r))   
        jproj = Int64(round(j*R/r))
        # copy likelihood from world to map
        p.observer.likelihood[Int64(i),Int64(j)] = p.observer.likelihood[iproj,jproj]
        p.observer.posterior[Int64(i),Int64(j)] = p.observer.posterior[iproj,jproj]
      end
    end
  end


end




 # Function to sample from normalized likelihood by rejection
 function sample_likelihood(p::Placozoan)

     n = 0
     while n < p.observer.nLparticles[]
       candidate =   2.0*p.observer.maxRange*(rand(2).-0.5) #rand(-p.observer.maxRange:p.observer.maxRange,2)
       r = sqrt(candidate[1]^2 + candidate[2]^2)
       if  r < p.observer.maxRange
         if rand() < p.observer.likelihood[Int64.(round.(candidate[1])), Int64.(round.(candidate[2]))]
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
function stalk(predator::Placozoan, prey::Placozoan, min_Δ::Float64)

# predator location in polar coordinates
predatorRange  = sqrt(predator.x[]^2 + predator.y[]^2)  # predator distance from origin

Δ = predatorRange - prey.radius - predator.radius  # distance between edges

driftDirection = sign(Δ - min_Δ)

if driftDirection > 0.0  # predator is approaching

  # model prey spinning by rotating predator trajectory
  predator.x[] =  cos(prey.ω[])*predator.x[] + -sin(prey.ω[])*predator.y[]
  predator.y[] =  sin(prey.ω[])*predator.x[] + cos(prey.ω[])*predator.y[]

end

#speed = predator.speed[]*(50.0/abs(Δ))^0.3
speed = predator.speed[]*(100.0/(50.0+Δ))^(2.0)

# drift towards prey, hold at min_Δ
predator.x[] -= driftDirection*speed*predator.x[]/predatorRange  # x = x - speed*sin(bearing) 
predator.y[] -= driftDirection*speed*predator.y[]/predatorRange  # y = y - speed*sin(bearing) 

# Gaussian diffusion
predator.x[] +=  0.75*randn()[] # x = x  + noise
predator.y[] +=  0.75*randn()[] # y = y  + noise 


end


function initialize_particles(p::Placozoan)

  p.observer.Pparticle[1:p.observer.nPparticles[],:] = samplePrior(p.observer.nPparticles[], p)


end

# function bayesParticleUpdate(p::Placozoan)

#   δ2 = p.observer.collisionRadius[]^2
#   #diffuseCoef = 2.0   # posterior particle diffusion rate (SD of Gaussian per step) at edge of mat
#   # NB diffusion coef here should match diffusion coef in sequential Bayes upsdate (bayesArrayUpdate())
#   nSpawn = 4  # average number of new posterior particles per collision
#   nCollision = 0
#   nCollider = 0
#   collision = fill(0, p.observer.nPparticles[])
#   collider = fill(0, p.observer.nPparticles[])


#   # diffuse (random gaussian jitter) posterior particles
#   # with boundary at edge of mat. Include vigilance bias by (rarely) replacing
#   # a posterior particle by a sample from the prior
#   @inbounds for i in 1:p.observer.nPparticles[]

#     if rand()[]<p.observer.posteriorDeathRate[] # with probability = posterior death rate
#       p.observer.Pparticle[i,:] = samplePrior(1,p)  # replace posteior particle with sample from prior
#     else   # otherwise particle Brownian motion
#       on_mat = false
#       while !on_mat

#         # compute a candidate location for each posterior particle to jump to 
#         candidate = p.observer.Pparticle[i,:] 
#         r = sqrt(candidate[1]^2 + candidate[2]^2)
#         candidate += p.observer.diffuseCoef[]*randn(2)*r/p.observer.maxRange # uniform diffusion in epithelium

#         # if the candidate location is within the boundary of the internal map
#         # then move the particle to that location. Otherwise try again.
#         r = sqrt(candidate[1]^2 + candidate[2]^2)
#         if (r > p.radius) & (r < p.observer.maxRange) # if candidate location is on the mat
#           p.observer.Pparticle[i,:] = candidate       # move particle to candidate location
#           on_mat = true                               # and we're done. Otherwise try again.
#         end
#       end  
#     end
#   end

#     # # On each update a fixed number (proportion) of posterior particles 
#     # # die at random and are reincarnated as (replaced by)
#     # #  a random sample from the initial prior.
#     # # (stops posterior particles prematurely condensing into local clouds, 
#     # #  maintains 360 deg attention; biophysically interpreted as equilibrium 
#     # #  between production and decay of posterior particles)
#     # #nscatter = Int(round(p.observer.priorDensity[]*p.observer.nPparticles[]))
#     # iscatter = rand(1:p.observer.nPparticles[], p.observer.posteriorDeaths[] )
#     # p.observer.Pparticle[iscatter, :] = samplePrior(p.observer.posteriorDeaths[], p)
#     # #p.observer.Pparticle_step[iscatter,:] .= 0.0

#   # list Pparticles that have collided with Lparticles
#      nL = p.observer.nLparticles[]
#      L = p.observer.Lparticle[:,:]
#      @inbounds for i = 1:p.observer.nPparticles[]  # find collisions between posterior and likelihood particles
#       #J = shuffle(1:nL)
#       @inbounds for j in 1:nL
#           if ((p.observer.Pparticle[i, 1] - L[j, 1])^2 +
#                 (p.observer.Pparticle[i, 2] - L[j, 2])^2) < δ2  #  collision
#           # if (abs((p.observer.Pparticle[i, 1] - L[j, 1]))<p.observer.collisionRadius[]) &
#           #   (abs((p.observer.Pparticle[i, 2] - L[j, 2]))<p.observer.collisionRadius[])
#               nCollision = nCollision + 1
#               collision[nCollision] = i     # ith posterior particle has collided with a likelihood particle
#               L[j:(nL-1)] = L[(j+1):nL]     # remove the Lparticle from list of available colliders
#               nL = nL - 1
#             break
#           end
#         end
#       end

#    #  print(nCollision)

#   if nCollision > 0

#     # copy posterior particles that have collided with sensory particles
#     particle = p.observer.Pparticle[collision[1:nCollision], :] 

#     # each collision spawns a Poisson-distributed number of new particles
#     n_newparticles = rand(Poisson(nSpawn), nCollision)


#     @inbounds for i = 1:nCollision
#        @inbounds for j = 1:n_newparticles[i]   # replace randomly chosen posterior particles with offspring of collision
#           p.observer.Pparticle[rand(1:p.observer.nPparticles[]), :]  = particle[i,:]
#         end
#     end

#   end

#   # reflect likelihood, posterior, likelihood particles and posterior particles into the marginal zone
#   reflectParticles!(p)  

# end

function bayesParticleUpdate(p::Placozoan)

  δ2 = 1.5   # squared collision maxRange
  sδ = 1.0  # scatter maxRange
  nCollision = 0

  # collision[i] counts the number of Lparticles that the 
  # ith Pparticle has collided with
  collision = fill(0, p.observer.nPparticles[])
  # list Bparticles that have collided with Lparticles
  for i = 1:p.observer.nPparticles[]
    for j = 1:p.observer.nLparticles[]  # check for collisions with belief
      if (p.observer.Pparticle[i, 1] - p.observer.Lparticle[j, 1])^2 +
         (p.observer.Pparticle[i, 2] - p.observer.Lparticle[j, 2])^2 < δ2   # collision
        nCollision = nCollision + 1
        collision[i] += 1   
      end
    end
  end

  if nCollision > 0
    # each collision produces a Poisson-distributed number of new particles
    # this may produce an excess of Pparticles that will need to be winnowed
    # (that's why this array, which will hold the locations of new particles,
    #  is much larger than the number of new particles we will end up with)
    newParticle = fill(0.0, 20*p.observer.nPparticles[], 2)

  #######################################
  #   TBD: ensure all new particles are on the mat
  ################################################
    newparticlecount = 0 
    for i = 1:p.observer.nPparticles[]
      for j = 1:collision[i]  # for each collision of ith particle
        nnew = 1+rand(Poisson(4))  # number of particles created by collision
        for k = 1:nnew
          newparticlecount += 1
          R = Inf
          while (R > p.observer.maxRange) | (R < p.radius)  # ensure new particle is on the mat
              # create new particle near ith Pparticle
            newParticle[newparticlecount,:] = p.observer.Pparticle[i, :] + sδ * randn(2)
            R = sqrt(newParticle[newparticlecount,1]^2 + newParticle[newparticlecount,2]^2)
          end
        end
      end
    end

    # choose up to p.observer.nPparticles[] of these new particles at random
    k = shuffle(1:newparticlecount)  # random index into new particles
    if newparticlecount > p.observer.nPparticles[]
      newparticlecount = p.observer.nPparticles[]
    end
    # random selection of new particles (in random order)
    newParticle =   newParticle[k[1:newparticlecount],:]

    # replace randomly chosen existing particles with new particles
    for i = 1:newparticlecount
      p.observer.Pparticle[rand(1:p.observer.nPparticles[]), :]  = newParticle[i,:]
    end

    # diffusion
    for i = 1:p.observer.nPparticles[]
      R = Inf
      while (R > p.observer.maxRange) | (R < p.radius) 
        p.observer.Pparticle[i,:] += 4.0*randn(2)
        R = sqrt(newParticle[newparticlecount,1]^2 + newParticle[newparticlecount,2]^2)
      end
    end

    # # kluge number of particles (normalize the discrete distribution)
    # # by random particle duplication
    # for i in 1:(p.observer.nPparticles[]-count)
    #   newBelief[count+i,:] = newBelief[rand(1:count),:]
    # end

    # p.observer.Pparticle[:] = newBelief[:]

    # draw 100S% of Bparticles from prior
    S = 0.001
    nscatter = Int(round(S * p.observer.nPparticles[]))
    iscatter = shuffle(1:p.observer.nPparticles[])[1:nscatter]
    p.observer.Pparticle[iscatter, :] = samplePrior(nscatter, p)
    p.observer.Pparticle_step[iscatter, :] = zeros(nscatter,2)
    # ϕ = 2 * π * rand(nscatter)
    # for i = 1:nscatter
    #   r = Inf
    #   while r > p.observer.maxRange[]
    #     r = p.observer.priormean[] + p.observer.priorsd[] * randn()[]
    #   end
    #   p.observer.Pparticle[iscatter[i], :] = [r * cos(ϕ[i]), r * sin(ϕ[i])]
    #   p.observer.Pparticle_step[iscatter[i], :] = [0.0, 0.0]
    # end

  end

end

# draw n samples from prior by rejection
function samplePrior(n, p::Placozoan)

  sample = zeros(n,2)
  count = 0
  top = maximum(p.observer.prior)

  while count < n
    
    # pick a random point on mat external to the placozoan
    r = p.radius +  (p.observer.maxRange - p.radius)*rand()
    ϕ = 2.0*π*rand()
    x = r*cos(ϕ)
    y = r*sin(ϕ)
   
    # accept candidate with probability proportional to prior at that point
    if top*rand() < p.observer.prior[Int(round(x)), Int(round(y))]
      count = count + 1
      sample[count,:] = [x,y]
    end
  
  end

  return (sample)

end


function initialize_prior(p::Placozoan)
# prior probability increases linearly with distance from origin
# (assume prey is more likely to appear further away)
  priorsum = 0.0
  priorMean = p.radius + 100.0        # expect predator to appear when it gets within 100um
  priorSigma2 = 100.0^2     # variance
  p.observer.prior .= 0.0
  for i = -p.observer.maxRange:p.observer.maxRange
    for j = -p.observer.maxRange:p.observer.maxRange
      d = sqrt(i^2 + j^2)
      if ( (d>p.radius) & (d<p.observer.maxRange) )
         p.observer.prior[i, j] = exp(-(d-priorMean)^2.0/priorSigma2)
         priorsum += p.observer.prior[i, j]
      end
    end
  end
   p.observer.prior[:,:] ./= priorsum
   p.observer.posterior[:,:] = p.observer.prior[:,:]  
end


function bayesArrayUpdate(p::Placozoan)

 
  # posteriorSum = 0.0
  # for i in -p.observer.maxRange:p.observer.maxRange
  #   for j in -p.observer.maxRange:p.observer.maxRange
  #       # d = sqrt(i^2 + j^2)
  #       # if (d>p.radius) & (d<p.observer.maxRange)
  #         # posterior is dynamic prior
  #         p.observer.posterior[i,j] *= p.observer.likelihood[i,j]
  #         posteriorSum += p.observer.posterior[i,j]
  #       # end
  #   end
  # end

  posteriorSum = sum(p.observer.posterior)
  p.observer.posterior .*= p.observer.likelihood/posteriorSum
 

   # diffuse and mix with initial prior
   # NB diffusion coef here should match diffusion coef in particle filter
   #density = p.observer.posteriorDeaths[]/p.observer.nPparticles[]
   #diffuseCoef = 2.0

   ###############################################################
   p.observer.posterior[:,:]  = (1.0-p.observer.posteriorDeathRate[])*
      imfilter(p.observer.posterior, Kernel.gaussian(diffuseCoef)) + p.observer.posteriorDeathRate[].*p.observer.prior[:,:]
      ###############################################################
      # p.observer.posterior[:,:]  = (1.0-density)*
      # imfilter(p.observer.posterior, Kernel.gaussian(diffuseCoef))./posteriorSum + density.*p.observer.prior[:,:]

  # renormalize over mat (ie compensate for probability mass that has leaked out of the observable world,
  # corresponding to not allowing particles to diffuse off the mat)
  # posteriorSum = 0.0
  # @inbounds for i in -p.observer.maxRange:p.observer.maxRange
  #   @inbounds for j in -p.observer.maxRange:p.observer.maxRange
  #       d = sqrt(i^2 + j^2)
  #       if (d>p.radius) & (d<p.observer.maxRange)
  #         posteriorSum += p.observer.posterior[i,j]
  #       else
  #         p.observer.posterior[i,j] = 0.0
  #       end
  #   end
  # end
  # p.observer.posterior[:,:] ./= posteriorSum

  posteriorSum = 0.0
for i in -p.observer.maxRange:p.observer.maxRange
for j in -p.observer.maxRange:p.observer.maxRange
        d = sqrt(i^2 + j^2)
        if (d>p.radius) & (d<p.observer.maxRange)
          posteriorSum += p.observer.posterior[i,j]
        else
          p.observer.posterior[i,j] = 0.0
        end
    end
  end
  p.observer.posterior[:,:] ./= posteriorSum

  # map likelihood and posterior density into placozoan marginal zone
  reflectArrays!(p)

end

# Utility functions

# make offset array entries radially symmetric around (0,0)
# i.e. make values depend only on distance from origin
function radialSmooth(X::OffsetArray, r::UnitRange{Int64})

  d = abs(X.offsets[1])  # array index offset 
  n = length(r)
  s = zeros(n)
  # compute average at each radius
  for i = 1:n
    count = 0   # count number of elements at radius r[i]
    for j in -d:d
      for k in -d:d
        if abs(sqrt(j^2 + k^2) - r[i])<1
          count = count + 1
          #s[i] = ((count-1)*s[i] + X[j,k])/count  # recursive mean
          s[i] += X[j,k]
        end
      end
    end
    s[i] = s[i]/count
  end

  # overwrite X with radial smoothed values
  for i = 1:n
    for j in -d:d
      for k in -d:d
        if abs(sqrt(j^2 + k^2) - r[i])<1
          X[j,k] = s[i]
        end
      end
    end
  end

end


# summarize particle distributions
function particleStats(prey::Placozoan, predator::Placozoan)

  bearing2predator = atan(predator.y[], predator.x[])    # bearing to centre of predator
  
  x_edge = predator.x[] - predator.radius*cos(bearing2predator) # x-coord of closest edge point
  y_edge = predator.y[] - predator.radius*sin(bearing2predator)

  bearing = atan(y_edge, x_edge)   # bearing to closest edge of predator
 

  # index active particles
  N = prey.observer.nPparticles[]

  # sorted distance from edge of prey to posterior particles
  D = sort(sqrt.(sum(prey.observer.Pparticle[1:N,:].^2, dims=2)), dims=1).-prey.radius

  # quantiles of particle distance to predator
  QN = [D[Int(round(q*N))] for q in  [0.01 0.05 0.25 0.50]]

  # estimated probability of predator within specfied range
  range = [40 45 50]
  NR = [sum(x->x<range[i], D)/N for i in 1:length(range)]

  # bearing error for each particle
  θ = atan.(prey.observer.Pparticle[1:N,2],prey.observer.Pparticle[1:N,1]) .- bearing

  # unwrap
  for i in 1:N
    if θ[i] > π
      θ[i] = θ[i] - 2π
    elseif    θ[i] < -π
      θ[i] = θ[i] + 2π
    end
  end

  # quantiles of particle bearing angle from predator
  # 1%, 5% amd 50% credibility intervals + median
  θ = sort(θ)
  Qθ = [ θ[Int(round(q*N))]*180/π for q in  [0.005 0.025 0.25 0.5 0.75 0.975 0.995] ]
  iQmed = 4  # index of median angle

  # M-cell threat estimate
  prey.mcell.x[] = prey.mcell.d*cos(bearing2predator + Qθ[iQmed]*π/180.0)
  prey.mcell.y[] = prey.mcell.d*sin(bearing2predator + Qθ[iQmed]*π/180.0)
  p = 0   # initialize particle-in-threat-zone count
  for i = 1:N
    if sqrt( (prey.mcell.x[]-prey.observer.Bparticle[i,1])^2 + 
             (prey.mcell.y[]-prey.observer.Bparticle[i,2])^2) <= prey.mcell.r
      p = p + 1
    end
  end
  p = p/N  # particle count to probability estimate

  
  # return proportion of particles in range, quantiles of particle range, 
  #   quantiles of particle angle and probability of predator in M-cell posterior field
  # ie M-cell's belief that there is a predator in its patch
  return (NR, QN, Qθ, p)

end

# summarize Bayesian observer distributions
function observerStats(prey::Placozoan, predator::Placozoan)

  bearing2predator = atan(predator.y[], predator.x[])    # bearing to centre of predator
  
  x_edge = predator.x[] - predator.radius*cos(bearing2predator) # x-coord of closest edge point
  y_edge = predator.y[] - predator.radius*sin(bearing2predator)

  bearing = atan(y_edge, x_edge)*180.0/π  # bearing to closest edge of predator

  #  radial cumulative distribution of posterior probability (ie integrate over direction)
  RCDF = fill(0.0, prey.observer.maxRange+1)
  # angular cumulative distn in 400 bins (ie using 400 instead of 360 "degrees")
  ACDF = fill(0.0, 361)
  # compute the radial density
  for i in -prey.observer.maxRange:prey.observer.maxRange
  for j in -prey.observer.maxRange:prey.observer.maxRange
      
    d = sqrt(i^2 + j^2)
    Θ = atan(j,i)*180/π - bearing  
    if Θ > 180.0
      Θ = Θ - 360.0
    elseif    Θ < -180.0
      Θ = Θ + 360.0
    end
    
    if (d< prey.observer.maxRange) & (d>=prey.radius)
      d = 1+Int64(round(d))     # quantized distance to world location
      RCDF[d] = RCDF[d] + prey.observer.posterior[i,j]
      Θ = 181 + Int64(round(Θ))      # quantized angle
      ACDF[Θ] =  ACDF[Θ] + prey.observer.posterior[i,j]
    end
    end
  end

  # cumulative range distribution 
  RCDF = cumsum(RCDF)
  
  # probability that predator is within specified range(s)
  # range values here should match those in particleStats()
  PR = RCDF[Int64(round(prey.radius)).+[40, 45, 50]]

  # quantiles of posterior density of distance to predator
  # quantile values here should match those in particleStats()
  QP = [ minimum(findall(x->x>=q, RCDF))-prey.radius for q in  [0.01 0.05 0.25 0.5] ]

  # cumulative angle distribution (clockwise re -y direction)
  ACDF = cumsum(ACDF)

  # angle quantiles (for credibility intervals)
  # 1% = 0.5% each end = 2/400 etc, note indexing from 1 not 0
  QΘ = [ minimum(findall(x->x>=q, ACDF)) for q in  [0.01 0.05 0.25 0.5 .75 .95 .99] ] .-180.0

  # update mcell coords
  prey.mcell.x[] = prey.mcell.d*cos(bearing2predator)
  prey.mcell.y[] = prey.mcell.d*sin(bearing2predator)
  # probability of predator in M-cell RF
  MP = posteriorInMcellRF(prey)

  (PR, QP, QΘ, MP)

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

# Entropy of empirical pdf in bits
# NB assumes sum(pdf)=1
function entropy(pdf)
  S = 0.0
  for p in pdf
    if p>1.0e-14  # prevent NaN error
      S = S - p*log2(p)
    end
  end
  S
end


#  Record the distance Δ between edges of predator and prey on each time step
function getRange(prey::Placozoan, predator::Placozoan, i)
    prey.observer.Δ[i] = sqrt(predator.x[]^2 + predator.y[]^2)-prey.radius-predator.radius
end

# # KLD(bayesian||particle)
# # Kullback-Liebler divergence
# # information gained by using particle distribution instead of bayesian posterior
# # i.e.  overconfidence in the particle distribution
# function KLD!(I::Observer, frame::Int64)

#   # KLD of particle estimate
#   S = 0.0
#   n = 0
#   # at each step the posterior includes a 
#   outlier_threshold = 1.0e-6
#   @inbounds for k in 1:I.nPparticles[]
#     i = Int64(round(I.Pparticle[k,1]))
#     j = Int64(round(I.Pparticle[k,2]))
#     if I.posterior[i,j] > outlier_threshold
#       S = S + log2(I.posterior[i,j])
#       n = n + 1
#     end
#     #end
#   end
#   I.KLD[frame] = S/n +  log2(n)

#   # KLD of random uniform sample
#   S0 = 0.0
#   nSamples = 0
#   n = 0
#   while nSamples < I.nPparticles[]
#     i = rand(-I.maxRange:I.maxRange,1)[]
#     j = rand(-I.maxRange:I.maxRange,1)[]
#     d = sqrt(i^2+j^2)
#     if (d>=I.minRange) & (d<=I.maxRange) # exclude particles not in the observable world
#       nSamples = nSamples + 1
#       if I.posterior[i,j] > outlier_threshold
#        #S0 = S0 + I.posterior[i,j]*log2(I.posterior[i,j] )
#        S0 = S0 + log2(I.posterior[i,j] )
#        n = n + 1
       
#       end
#     # end
#     end
#   end
#   I.KLD0[frame] = S0/n +  log2(n)


#    # KLD of sample from posterior
#    SI = 0.0
#    s = sample(I.posterior, I.nPparticles[])
#    n = 0
#    @inbounds for i in 1:I.nPparticles[]
#        #SI = SI + I.posterior[s[i,1],s[i,2]]*log2(I.posterior[s[i,1],s[i,2]])
#        if I.posterior[s[i,1],s[i,2]] > outlier_threshold
#         SI = SI + log2(I.posterior[s[i,1],s[i,2]])
#         n = n + 1
#        end
#    end

#    I.KLDI[frame] = SI/n +  log2(n)

#  end

# KLD(particle||bayesian)
function KLD!(I::Observer, frame::Int64)
  # computes Kullback-Liebler divergence, saves in KLD field of observer.
  # KLD is the information lost if the particle distribution is used instead of the 
  # "true" posterior. It is computed for the placozoan model (KLD), a null model in 
  # which particles are drawn uniformly at random over the mat (KLD0) and an "ideal"
  # particle distribution drawn from the true posterior (KLDI).
  # THIS VERSION of KLD compares the true posterior at sample points (particle locations)
  # to the particle distribution i.e. a uniform distribution over those points.

  # KLD of particle estimate
  S = 0.0
  n = 0
  outlier_threshold = 1.0e-14
  psum = 0.0
  @inbounds for k in 1:I.nPparticles[]
    i = Int64(round(I.Pparticle[k,1]))
    j = Int64(round(I.Pparticle[k,2]))
    try  # rarely [i,j] is out of bounds, but we dont care
       psum += I.posterior[i,j]
    catch
    end
  end

  @inbounds for k in 1:I.nPparticles[]
    i = Int64(round(I.Pparticle[k,1]))
    j = Int64(round(I.Pparticle[k,2]))
    #if (i^2 + j^2)<I.maxRange^2 # exclude particles not in the observable world
    try
      if I.posterior[i,j] > outlier_threshold
        # S = S + I.posterior[i,j]*log2(I.posterior[i,j] )
        S = S + I.posterior[i,j]*log2(I.posterior[i,j]/psum) 
      end
    catch
    end
  end
  I.KLD[frame] = S/psum + log2(I.nPparticles[])

  # KLD of random uniform sample
  S0 = 0.0
  nSamples = 0
  p0sum = 0.0
  while nSamples < I.nPparticles[]
    i = rand(-I.maxRange:I.maxRange,1)[]
    j = rand(-I.maxRange:I.maxRange,1)[]
    d = sqrt(i^2+j^2)
    if (d>=I.minRange) & (d<=I.maxRange) # exclude particles not in the observable world
      nSamples = nSamples + 1
      if I.posterior[i,j] > outlier_threshold
       p0sum = p0sum + I.posterior[i,j]
      end
    end
  end
  nSamples = 0
  while nSamples < I.nPparticles[]
    i = rand(-I.maxRange:I.maxRange,1)[]
    j = rand(-I.maxRange:I.maxRange,1)[]
    d = sqrt(i^2+j^2)
    if (d>=I.minRange) & (d<=I.maxRange) # exclude particles not in the observable world
      nSamples = nSamples + 1
      if I.posterior[i,j] > outlier_threshold
       #S0 = S0 + I.posterior[i,j]*log2(I.posterior[i,j] )
       S0 = S0 + I.posterior[i,j]*log2(I.posterior[i,j]/p0sum)
      end
    end
  end
  I.KLD0[frame] = S0/p0sum +  log2(I.nPparticles[])


   # KLD of sample from posterior
   SI = 0.0
   s = sample(I.posterior, I.nPparticles[])
   n = 0
   pIsum = 0.0
   @inbounds for i in 1:I.nPparticles[]
       #SI = SI + I.posterior[s[i,1],s[i,2]]*log2(I.posterior[s[i,1],s[i,2]])
       if I.posterior[s[i,1],s[i,2]] > outlier_threshold
        pIsum = pIsum + I.posterior[s[i,1],s[i,2]]
       end
   end

   @inbounds for i in 1:I.nPparticles[]
    #SI = SI + I.posterior[s[i,1],s[i,2]]*log2(I.posterior[s[i,1],s[i,2]])
    if I.posterior[s[i,1],s[i,2]] > outlier_threshold
     SI = SI + I.posterior[s[i,1],s[i,2]]*log2(I.posterior[s[i,1],s[i,2]]/pIsum)
    end
end
   I.KLDI[frame] = SI/pIsum +   log2(I.nPparticles[])

 end


function sample!(s::Array{Int64, 2}, D::AbstractArray)
  # draw sample s (nx2) from 2D empirical distribution D by rejection
  # sum(D)==1.
  # samples are returned as Int64 nx2 indices of D
  # Uses rectangular uniform proposal distribution whose 1/2-width shrinks
  # to 3x the standard deviation of the posterior as it converges. 

 (peak, ipeak) = findmax(D)   # maximum probability and its location
  σ = 2^((entropy(D)-4.0942)/2)  # standard deviation of target distribution

  N = size(s,1)      # required sample size
  i = 0              # sample size counter
  X = axes(D,1)
  Y = axes(D,2)

  # sample region is 3 sd each side of peak
  Δ = Int64(round(4*σ))
  x0 = max( minimum(X),   ipeak[1] - Δ)
  x1 = min( maximum(X), ipeak[1] + Δ)
  y0 = max( minimum(Y),   ipeak[2] - Δ)
  y1 = min( maximum(Y), ipeak[2] + Δ) 

  while i<N
    x = rand(x0:x1,1)[]    # uniform random point
    y = rand(y0:y1,1)[]
    if peak*rand()[] < D[x,y]  # accept/reject
      i = i + 1
      s[i,1] = x
      s[i,2] = y
    end
  end

  s
end

function sample(D::AbstractArray, N::Int64)
  # draw sample s of size n from 2D empirical distribution D by rejection]
  # sum(D)==1.
  # samples are returned as Int64 nx2 indices of D
  # samples are returned as Int64 nx2 indices of D
  # Uses rectangular uniform proposal distribution whose 1/2-width shrinks
  # to 3x the standard deviation of the posterior as it converges. 

  s = fill(0, N, 2)   # for samples

  (peak, ipeak) = findmax(D)  # maximum probability and its location
  σ = 2^((entropy(D)-4.0942)/2)  # standard deviation of target distribution


  i = 0              # sample size counter
  X = axes(D,1)
  Y = axes(D,2)

  # sample region is 3 sd each side of peak
  Δ = Int64(round(3*σ))
  x0 = max(minimum(X), ipeak[1] - Δ)
  x1 = min(maximum(X), ipeak[1] + Δ)
  y0 = max(minimum(Y), ipeak[2] - Δ)
  y1 = min(maximum(Y), ipeak[2] + Δ) 

  while i<N
    x = rand(x0:x1,1)[]    # uniform random point
    y = rand(y0:y1,1)[]
    if peak*rand()[] < D[x,y]  # accept/reject
      i = i + 1
      s[i,1] = x
      s[i,2] = y
    end
  end

  s
end

# posterior probability of predator in M-cell RF
function posteriorInMcellRF(p::Placozoan)
   
  # get m-cell vertices
  mc_pts = decompose(Point2f0, Circle(Point2f0(p.mcell.x[],p.mcell.y[]), p.mcell.r))
  mc_pts[end] = mc_pts[1]  # close polygon
  rf_pts = copy(mc_pts)  

  # project to get RF vertices, 
  # and record the min and maxy x and y coords of the RF
  xMin = yMin = Inf
  xMax = yMax = -Inf
  for j in 1:64
    pt = mc_pts[j]
    Ω = atan(pt[2], pt[1])
    r = sqrt(pt[1]^2 + pt[2]^2)
    r1 = (p.radius - r)*(p.observer.maxRange - p.radius)/p.marginwidth + p.radius
    x = r1*cos(Ω)
    y = r1*sin(Ω)
    rf_pts[j] = Point2f0(x,y)
    if x<xMin xMin = x end
    if x>xMax xMax = x end
    if y<yMin yMin = y end
    if y>yMax yMax = y end
  end

  # close polygon (required by inpolygon() )
  rf_pts[end] = rf_pts[1]

  # integrate probability in RF
  j0 = Int64(floor(yMin))
  j1 = Int64(ceil(yMax))
  Pr = 0.0
  for i in Int64(floor(xMin)):Int64(ceil(xMax))
 for j in j0:j1
      if inpolygon(Point2f0(i,j), rf_pts)==1
        #Pr = Pr + p.observer.posterior[i,j]
        Pr += p.observer.posterior[i,j] 
      end
    end
  end
 
  Pr
end

