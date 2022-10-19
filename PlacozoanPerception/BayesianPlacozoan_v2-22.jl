# BayesianPlacozoan module
# types and functions &c used by PlacozoanStalker.jl

using GLMakie
using Colors
using OffsetArrays
using Distributions
using Random
using ImageFiltering
using CSV
using DataFrames
using PolygonOps  

mcell_radius = 2.5

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



# Particle sizes  
const size_likelihood = 2.0
#size_prior = 4
const size_posterior = 2.5

const size_observation = 1.5
#size_prediction = 2
const size_belief = 1.5

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
  position::Vector{Point2{Float64}}
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
  Lparticle::Vector{Point2{Float64}}       # likelihood particles 
  Pparticle::Vector{Point2{Float64}}    # posterior particles
  Sparticle::Vector{Point2{Float64}}    # sensory  particles (reflected likelihood particles)
  Bparticle::Vector{Point2{Float64}}    # belief particles (reflected posterior particles)
  Pparticle_step::Vector{Point2{Float64}}  # particle prediction steps
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
               fill(Point2(0.0, 0.0), nPparticles),
               fill(Point2(0.0, 0.0), nPparticles),
               fill(Point2(0.0, 0.0), nPparticles),
               fill(Point2(0.0, 0.0), nPparticles),
               fill(Point2(0.0, 0.0), nPparticles),
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
  z2 = fill(Point2(0.0, 0.0), 1)
  zOff = zeros(-1:1, -1:1)
  Observer(1, 1, 1, 1.0, zOff, zOff, zOff, [1], [1], z2, z2, z2, z2, z2, 
          [0.0],[0.0], [0.0], 0, z1, z1, z1, z1, z1)
end



# Electroreceptor definition
struct Ereceptor
  N::Int64
  size::Float64  # symbol size for drawing receptors
  position::Vector{Point2{Float64}}
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
   position = [Point2(placozoanradius.*(cos(2π*i/N)), placozoanradius.*(sin(2π*i/N))) for i in 1:N ]

   # 1d vector containing N offset arrays; ith will contain RF for ith receptor
   Lhd = Array{OffsetArray,1}(undef,N)
   for i in 1:N
     Lhd[i] = zeros(-worldradius:worldradius,-worldradius:worldradius)
   end

   return Ereceptor(N, receptorSize, position, zeros(N), Lhd, openColor, closedColor)
end

# dummy electroreceptor constructor (for constructing placozoan without electroreceptors)
function Ereceptor()

  return Ereceptor(0, 0, [Point2(0.0, 0.0)], zeros(1),
                 Array{OffsetArray,1}(undef,1), RGB(0,0,0), RGB(0,0,0))
end


# Crystal cell receptor definition
struct CrystalCell
  N::Int64
  size::Float64  # symbol size for drawing receptors
  position::Array{Point2{Float64},1}  
  state::Array{Float64,1} # 0/1 for receptor in closed/open state
  lineOfSight::Array{Float64,1}
  pOpenV::Array{OffsetArray,1}
  lightColor::RGB
  darkColor::RGB
end

# crystal cell constructor
function CrystalCell(worldradius::Int64, xtalradius::Float64,
                   N::Int64, crystalSize::Float64,
                   lightColor::RGB, darkColor::RGB)
   if floor(N/4)!=N/4  ###???floor?? code for lowest?
     error("Number of receptors must be a multiple of 4")
   end

   # N receptors equally spaced in a ring at radius of (gut?/Inside margin?)
   lineOfSight = [2π*i/N for i in 0:(N-1)]  # each xtal faces radially outwards
   position = xtalradius*Point2.(cos.(lineOfSight), sin.(lineOfSight))
   # lineOfSight = [atan(y[i],x[i]) for i in 1:N]
   # 1d vector containing N offset arrays; ith will contain RF for ith receptor
   Lhd = Array{OffsetArray,1}(undef,N)
   for i in 1:N
     Lhd[i] = zeros(-worldradius:worldradius,-worldradius:worldradius)
   end

   return CrystalCell(N, crystalSize, position, zeros(N), lineOfSight,
                    Lhd, lightColor, darkColor)
end

# crystal cell dummy constructor
function CrystalCell()

  return CrystalCell(0, 0, [Point2(0.0,0.0)], zeros(1), [0],
                 Array{OffsetArray,1}(undef,1), RGB(0,0,0), RGB(0,0,0))
end

# Placozoan definition
struct Placozoan
  radius::Float64
  marginwidth::Float64
  gutradius::Float64
  celldiam::Float64
  position::Vector{Point2{Float64}}  # mutable position[][i]
  # field[i] is pre-computed bio-electric field strength
  #   at distance i μm from edge of placozoan
  field::Array{Float64,1}
  potential::Array{Float64,1}  # in μV
  fieldrange::Int64   # number of elements in field (= max maxRange in μm)
  receptor::Ereceptor  # electroreceptor array
  photoreceptor::CrystalCell
  observer::Observer  # Bayesian particle filter
  mcell::Mcell        # "Mauthner" neuron
  pcell::Mcell        # particle M-cell
  speed::Array{Float64,1}  # speed parameter (actual speed depends on state)
  brown::Array{Float64,1}  # Brownian motion parameter (diffusion constant)
  ω::Array{Float64,1}   # omega, angular velocity of self-rotation
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
      [Point2(0.0, 0.0)],
      zeros(fieldrange),
      zeros(fieldrange),
      fieldrange,
      receptor,
      crystalcell,
      observer,
      Mcell(radius-7.14, mcell_radius, [Point2(0.0, 0.0)] ),
      Mcell(radius-5.5, mcell_radius, [Point2(0.0, 0.0)] ),
      [0.0],  
      [diffuseCoef],
      [0.0],
      [0.0, 0.0],
      bodycolor,
      gutcolor,
      edgecolor,
    )

end # Placozoan constructor

# placozoan constructor with field, initially not at the origin,
# it moves but but no receptors or observer
# i.e. a predator in the initial scenario

function Placozoan(radius::Int64, margin::Int64, fieldrange::Int64, position::Point2{Float64}, speed::Float64, brownian::Float64,
                  bodycolor::RGBA, gutcolor::RGBA, edgecolor::RGB)

   return Placozoan(radius, margin, radius-margin, 12.0,  [position],
     zeros(fieldrange), zeros(fieldrange), fieldrange,
     Ereceptor(), CrystalCell(), Observer(), 
     Mcell(0.0, 0.0, [Point2(0.0, 0.0)]), Mcell(0.0, 0.0, [Point2(0.0, 0.0)]),
     [speed], [brownian], [0.0], [0.0, 0.0],
     bodycolor, gutcolor, edgecolor )

end

# function initializeObserver(p::Placozoan, nLparticles::Int64, nPparticles::Int64,
#   priorDensity::Float64)

#    p.observer.nLparticles[]  = nLparticles
#    p.observer.nPparticles[]  = nPparticles
#    p.observer.priorDensity[] = priorDensity

# # +   p.observer.Lparticle = zeros(nLparticles,2)
# #    p.observer.Pparticle = zeros(nPparticles,2)
# #    p.observer.Pparticle_step = zeros(nPparticles,2)


#    likelihood(p)           # initialize likelihood given initial receptor states
#    sample_likelihood(p)    # sample from normalized likelihood
#    initialize_particles(p) # draw initial sample from prior
#    initialize_prior(p)     # initialize numerical Bayesian prior

# end


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
    # (because self must be still alive to do this computation :-)
    x = self.receptor.position[i][1]
    y = self.receptor.position[i][2]

    for j in -self.observer.maxRange:self.observer.maxRange
      for k in -self.observer.maxRange:self.observer.maxRange

        # likelihood at (j,k)
        L = sqrt(j^2+k^2) > self.radius ?
        Electroreceptor_pOpen(sqrt((x-j)^2 + (y-k)^2), other.potential) : 0.0
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
    x = self.photoreceptor.position[i][1]
    y = self.photoreceptor.position[i][2]    
    for j in -self.observer.maxRange:self.observer.maxRange
      for k in -self.observer.maxRange:self.observer.maxRange

        angleFromLineOfSight = atan(k-y,j-x) - 2π*(i-1)/self.photoreceptor.N

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

  # likelihood particles -> sensory particles
  R = distance(p.observer.Lparticle)
  r = ( p.radius .- p.marginwidth*(R.-p.radius)./(p.observer.maxRange-p.radius) )::Array{Float64,1}
  p.observer.Sparticle[:] = r.*p.observer.Lparticle[:]./R

  # posterior particles -> belief particles
  Rp = distance(p.observer.Pparticle)
  rp = (p.radius .- p.marginwidth*(Rp.-p.radius)./ (p.observer.maxRange-p.radius))::Array{Float64,1}
  p.observer.Bparticle[:] = rp.*p.observer.Pparticle[:]./Rp


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




 # sample from likelihood (by rejection)
 # nb maximum(likelihood)==1 
 function sample_likelihood(p::Placozoan)

  n = 1
  while n <= p.observer.nLparticles[]
    p.observer.Lparticle[n] = uniformSampleAnnulus(p.radius, Float64(p.observer.maxRange)) # random point on mat
    if rand() < p.observer.likelihood[gridPoint(p.observer.Lparticle[n])...]
      n = n + 1
    end
  end
 end


 # electroreceptor states as function of predator location
 function electroreception(prey::Placozoan, predator::Placozoan)
   
  for i = 1:length(prey.receptor.state)  

    range = max(distance(predator.position[] - prey.receptor.position[i]) - predator.radius, 0.0) 
    prey.receptor.state[i] = rand() < Electroreceptor_pOpen(range, predator.potential) ? 1 : 0

   end
 end


 # photoreceptor states as function of predator location (approach angle)
 function photoreception(prey::Placozoan, predator::Placozoan) 

  for j = 1:length(prey.photoreceptor.state)

     angle2predator = atan( predator.position[2] - prey.photoreceptor.position[j][2],
                               predator.position[1] - prey.photoreceptor.position[j][1])

     deviationFromLineOfSight = angle2predator - prey.photoreceptor.lineOfSight[j]

    #println(j, ", ", prey.photoreceptor.lineOfSight[j], ", ", angle2predator, ", ", deviationFromLineOfSight)

    #  if deviationFromLineOfSight > π
    #    deviationFromLineOfSight -= 2π
    #  end
    #  if deviationFromLineOfSight < -π
    #    deviationFromLineOfSight += 2π
    #  end

     prey.photoreceptor.state[j] = rand() < Photoreceptor_pOpen(deviationFromLineOfSight) ? 1 : 0

  end

end


# rotate placozoan location through angle dψ around origin
function orbit(dψ::Float64, placozoan::Placozoan)
  x = placozoan.position[][1]
  placozoan.position[] =  Point2( cos(dψ)*x + sin(dψ)*placozoan.position[][2], 
                                 -sin(dψ)*x + cos(dψ)*placozoan.position[][2])
  # C = [cos(dψ) sin(dψ); -sin(dψ) cos(dψ)]
end

# # rotate particles through angle dψ around origin
# function orbit(dψ::Float64, p::Array{Float64,2})
#   # p.x[] =  cos(dψ)*p.x[] + sin(dψ)*p.y[]
#   # p.y[] = -sin(dψ)*p.x[] + cos(dψ)*p.y[]
#   p = p*[cos(dψ) -sin(dψ); sin(dψ) cos(dψ)]
#    nothing
# end

# predator step
# returns Δ = distance between edges after step
function stalk(predator::Placozoan, prey::Placozoan, min_Δ::Float64)

  predatorRange = distance(predator.position)[]
  Δ = predatorRange - prey.radius - predator.radius  # distance between edges

  if Δ > 2.0*min_Δ  # predator not there yet

    # predator orbits prey (in prey frame, because prey is spinning in predator frame)
    # nb ω is omega not w
    orbit(prey.ω[], predator) 

  end

  # prey accelerates as it approaches
  #speed = predator.speed[]*(100.0/(50.0+Δ))

  # approach prey
  # nb predator.position[]/predatorRange  is a unit vector pointing from predator to prey
  # the predator will approach and hover at min_Δ
  predator.position[] -= sign(Δ - min_Δ)*predator.speed[]*predator.position[]/predatorRange 
 
  # diffusion/brownian motion: add Gaussian noise
  predator.position[] += predator.brown[]*Point2(randn(2))

   return(Δ)
end

# not needed anymore, samplePrior initializes Pparticles
# function initialize_p(p::Placozoan)

#   placozoan.observer.Pparticle[:] = samplePrior(placozoan)

# end


function bayesParticleUpdate(placozoan::Placozoan)

  #δ2 = 1.5   # squared collision maxRange
  #sδ = 1.0  # scatter maxRange
  nCollision = 0

  # collision[i] counts the number of Lparticles that the 
  # ith Pparticle has collided with
  collision = fill(0, placozoan.observer.nPparticles[])
  # list Bparticles that have collided with Lparticles
  # for i = 1:placozoan.observer.nPparticles[]
  #   for j = 1:placozoan.observer.nLparticles[]  # check for collisions with belief
  #     if distance(placozoan.observer.Pparticle[i], placozoan.observer.Lparticle[j]) < placozoan.observer.collisionRadius[]
  #       # jth likleihood particle has collided with ith posterior particle
  #       nCollision = nCollision + 1
  #       collision[i] += 1   
  #     end
  #   end
  # end

  i = 0
  for P in placozoan.observer.Pparticle
    i = i + 1
    # d = distance(P)
    # r = (d < placozoan.radius + 50.) ? placozoan.observer.collisionRadius[]*sqrt(distance(P)/placozoan.radius) : 1.2*placozoan.observer.collisionRadius[]
    for L in placozoan.observer.Lparticle # check for collisions with belief
      if distance(P,L) < placozoan.observer.collisionRadius[]*projectMap(P, placozoan)
        # jth likleihood particle has collided with ith posterior particle
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
    # we would get more than maxNewParticles only when Pparticles are so close to each other that 
    #   there is no point in looking at all of them
    maxNewParticles = 20*placozoan.observer.nPparticles[]  
    newParticle = fill(Point2(0.0, 0.0), maxNewParticles)

  #######################################
  #   TBD: ensure all new particles are on the mat
  ################################################
    newparticlecount = 0 
    OVERFLOW = false
    for i = 1:placozoan.observer.nPparticles[]
      for j = 1:collision[i]  # for each collision of ith particle
        nnew = 4 #1+rand(Poisson(2)) #*(distance(placozoan.observer.Pparticle[i])/placozoan.radius)^.5))  # number of posterior particles after collision 
        for k = 1:nnew
          newparticlecount += 1
          if newparticlecount > maxNewParticles
            # println("B")
            OVERFLOW = true
            newparticlecount = maxNewParticles
            break  # may have ot break multiple times, but don't care
          end
          R = Inf
          while (R > placozoan.observer.maxRange) | (R < placozoan.radius)  # ensure new particle is on the mat
              # create new particle near ith Pparticle
            newParticle[newparticlecount] = 
              placozoan.observer.Pparticle[i] + placozoan.observer.diffuseCoef[] * Point2(randn(2))
            R = distance(newParticle[newparticlecount])
          end
        end
        if OVERFLOW break end
      end
      if OVERFLOW break end
    end

    # choose up to p.observer.nPparticles[] of these new particles at random
    newParticle = newParticle[1:newparticlecount]
    shuffle!(newParticle)  # new particles in random order
    if newparticlecount > placozoan.observer.nPparticles[]
      newparticlecount = placozoan.observer.nPparticles[]
      newParticle =   newParticle[1:newparticlecount]
    end

    # replace randomly chosen existing particles with new particles
    shuffle!(placozoan.observer.Pparticle)
    for i = 1:newparticlecount
      placozoan.observer.Pparticle[i]  = newParticle[i]
    end

  end

  # diffusion
  for i = 1:placozoan.observer.nPparticles[]
    R = Inf
    while (R >placozoan.observer.maxRange) | (R < placozoan.radius) # retry if not on mat
      global candidate = placozoan.observer.Pparticle[i] + placozoan.observer.diffuseCoef[]*Point2(randn(2))
      R = distance(candidate)
    end
    placozoan.observer.Pparticle[i] = candidate
  end

    # draw 100S% of Bparticles from prior
  #S = 0.001
  NpriorParticles = Int(round(placozoan.observer.posteriorDeathRate[] * placozoan.observer.nPparticles[]))  # number of particles to draw & replace
  priorPoint = uniformSampleAnnulus(Float64(placozoan.observer.minRange), 
                    Float64(placozoan.observer.maxRange), NpriorParticles)# uniform random sample

  
  priorPoint = samplePrior(placozoan, NpriorParticles)   # random sameple from prior
  # println("e")
  shuffle!(placozoan.observer.Pparticle)                # randomize order of posterior points
  # println("f")
  for i in 1:NpriorParticles                             # replace randomly selected posterior points with samples from prior
    placozoan.observer.Pparticle[i] = priorPoint[i]
    placozoan.observer.Pparticle_step[i] = Point2(0.0, 0.0) # legacy - not used??
  end

  reflectParticles!(placozoan)

end

#  initialize posterior particles as random sample from prior (by rejection)
function initializePosteriorParticles(placozoan::Placozoan)

  maxPrior = maximum(placozoan.observer.prior)
  n = 1
  while n <= placozoan.observer.nPparticles[]
    placozoan.observer.Pparticle[n] = uniformSampleAnnulus(placozoan.radius, Float64(placozoan.observer.maxRange)) # random point on mat
    if maxPrior*rand() < placozoan.observer.prior[gridPoint(placozoan.observer.Pparticle[n])...]
      n = n + 1
    end
  end   
end

#  random sample of N particles from prior (by rejection)
function samplePrior(placozoan::Placozoan, N::Int64)

  maxPrior = maximum(placozoan.observer.prior)
  samplePoint = fill(Point2(0.0, 0.0), N)
  n = 1
  while n <= N
    samplePoint[n] = uniformSampleAnnulus(placozoan.radius, Float64(placozoan.observer.maxRange)) # random point on mat
    # if placozoan.observer.prior[gridPoint(placozoan.observer.Lparticle[n])...] < .001
    #   println(gridPoint(placozoan.observer.Lparticle[n]))
    # end
    if maxPrior*rand() < placozoan.observer.prior[gridPoint(samplePoint[n])...]
      n = n + 1
    end
  end   
  return(samplePoint)
end


# compute prior on grid, copy to initial posterior
function initializePosteriorPDF(p::Placozoan)

  priorsum = 0.0
  priorMean = p.radius + 175.0        # expect predator to appear when it gets within 100um
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

  posteriorSum = sum(p.observer.posterior)
  p.observer.posterior .*= p.observer.likelihood/posteriorSum
 

  # diffuse and mix with uniform prior 
  p.observer.posterior[:,:]  = (1.0-p.observer.posteriorDeathRate[])*
    imfilter(p.observer.posterior, Kernel.gaussian(p.observer.diffuseCoef[])) .+ 
              p.observer.posteriorDeathRate[]/p.observer.N

  # renormalize 
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

  # Θ = atan(predator.position[][2], predator.position[][1])    # bearing to centre of predator

  # point on the boundary of the predator that is closest to the prey
  pClosest = predator.position[] - predator.radius*(predator.position[]-prey.position[])/distance(predator.position[]-prey.position[])  

  # bearing angle to closest point
  bearing = atan(pClosest[2], pClosest[1])   
 

  # index active particles
  N = prey.observer.nPparticles[]

  # sorted distance from closest point to all posterior particles
  D = sort(distance(prey.observer.Pparticle[:])) .- prey.radius

  # quantiles of particle distance to predator
  QN = [D[Int(cround(q*N))] for q in  [0.01 0.05 0.25 0.50]]

  # estimated probability of predator within specfied range
  range = [40 45 50]
  NR = [sum(x->x<range[i], D)/N for i in 1:length(range)]

  
  #θ = atan.(prey.observer.Pparticle[:][2], prey.observer.Pparticle[:][1]) .- bearing

  # bearing error for each particle
  θ = fill(0.0, N)
  for i in 1:N
    θ[i] = atan(prey.observer.Pparticle[i][2], prey.observer.Pparticle[i][1]) .- bearing
    if θ[i] > π
      θ[i] = θ[i] - 2π
    elseif    θ[i] < -π
      θ[i] = θ[i] + 2π
    end
  end

  # quantiles of particle bearing angle from predator
  # 1%, 5% amd 50% credibility intervals + median
  θ = sort(θ)
  Qθ = [ θ[Int(cround(q*N))]*180/π for q in  [0.005 0.025 0.25 0.5 0.75 0.975 0.995] ]
  iQmed = 4  # index of median angle

  # M-cell threat estimate
  # nb M-cell is positioned dynamically on the predator bearing
  # so 1 M-cell simulates the closest of an array of M-cells
  prey.pcell.position[] = prey.pcell.d*Point2(cos(bearing + Qθ[iQmed]*π/180.0), sin(bearing + Qθ[iQmed]*π/180.0))
  p = 0   # initialize particle-in-threat-zone count
  for i = 1:N
    if distance(prey.pcell.position[], prey.observer.Bparticle[i]) <= prey.pcell.r
      p = p + 1
    end
  end
  p = p/N  # particle count to probability estimate

  
  # return proportion of particles in range, quantiles of particle range, 
  #   quantiles of particle angle and probability of predator in M-cell posterior field
  # ie M-cell's belief that there is a predator in its patch
  return (NR, QN, Qθ, p)

end

# Max posterior probability of placozoan location
function MAPlocation(placozoan::Placozoan)

  r = placozoan.radius
  D = placozoan.observer.maxRange
  p = findmax(OffsetArrays.no_offset_view(prey.observer.posterior).*[ ((i)^2+(j)^2)>(placozoan.radius^2) for i in -D:D, j in -D:D])
  
  # transform back to offset indices
  return(Tuple(p[2]) .- D .- 1)

end


# summarize Bayesian observer distributions
function observerStats(prey::Placozoan, predator::Placozoan)

  bearing2predator = atan(predator.position[][2], predator.position[][1])    # bearing to centre of predator
  
  x_edge = predator.position[][1] - predator.radius*cos(bearing2predator) # x-coord of closest edge point
  y_edge = predator.position[][2] - predator.radius*sin(bearing2predator)

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
      d = 1+Int64(cround(d))     # quantized distance to world location
      RCDF[d] = RCDF[d] + prey.observer.posterior[i,j]
      Θ = 181 + Int64(cround(Θ))      # quantized angle
      ACDF[Θ] =  ACDF[Θ] + prey.observer.posterior[i,j]
    end
    end
  end

  # cumulative range distribution 
  RCDF = cumsum(RCDF)
  
  # probability that predator is within specified range(s)
  # range values here should match those in particleStats()
  PR = RCDF[Int64(cround(prey.radius)).+[40, 45, 50]]

  # quantiles of posterior density of distance to predator
  # quantile values here should match those in particleStats()
  QP = [ minimum(findall(x->x>=q, RCDF))-prey.radius for q in  [0.01 0.05 0.25 0.5] ]

  # cumulative angle distribution (clockwise re -y direction)
  ACDF = cumsum(ACDF)

  # angle quantiles (for credibility intervals)
  # 1% = 0.5% each end = 2/400 etc, note indexing from 1 not 0
  QΘ = [ minimum(findall(x->x>=q, ACDF)) for q in  [0.01 0.05 0.25 0.5 .75 .95 .99] ] .-180.0

  # update mcell coords
  prey.mcell.position[] = prey.mcell.d * Point2( cos(bearing2predator), sin(bearing2predator) )
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

# KLD(particle||bayesian)
function DKLD!(placozoan::Placozoan, frame::Int64)
  # Kullback-Liebler divergence from posterior distribution at points drawn from itself (posterior sample distribution)
  #   to posterior distribution at posterior particle locations, and KLD from posterior 
  #   sample distribution to posterior distribution at prior.

  # distribution of posterior at random sample from prior.
  priorSample = samplePrior(placozoan, placozoan.observer.nPparticles[])  # random sample from prior
  posteriorPDF_priorSample = fill(0.0,placozoan.observer.nPparticles[])
  for i in 1:placozoan.observer.nPparticles[]
    posteriorPDF_priorSample[i] = placozoan.observer.posterior[gridPoint(priorSample[i])...] # posterior density at sample points
  end
  posteriorPDF_priorSample /= sum(posteriorPDF_priorSample) # normalize

  # distribution of posterior at random sample from posterior.
  postSample = samplePosterior(placozoan)  # random sample from prior
  posteriorPDF_posteriorSample = fill(0.0,placozoan.observer.nPparticles[])
  for i in 1:placozoan.observer.nPparticles[]
    posteriorPDF_posteriorSample[i] = placozoan.observer.posterior[gridPoint(postSample[i])...] # posterior density at sample points
  end
  posteriorPDF_posteriorSample /= sum(posteriorPDF_posteriorSample) # normalize
  

    # distribution of posterior at posterior particle locations
    priorSample = samplePrior(placozoan, placozoan.observer.nPparticles[])  # random sample from prior
    posteriorPDF_posteriorParticles = fill(0.0,placozoan.observer.nPparticles[])
    for i in 1:placozoan.observer.nPparticles[]
      posteriorPDF_posteriorParticles[i] = placozoan.observer.posterior[gridPoint(prey.observer.posterior[i])...] # posterior density at sample points
    end
    posteriorPDF_posteriorParticles  # normalize

end



# KLD(particle||bayesian)
function KLD!(placozoan::Placozoan, frame::Int64)
  # (1) Kullback-Liebler divergence from posterior particle distribution to posterior distribution
  # (2) ... from posterior sample to posterior
  # (3) ... from prior sample to posterior



  # KLD of particle estimate
  S = 0.0
  n = 0
  outlier_threshold = 1.0e-14
  psum = 0.0

  # sum of posterior densitoes at sample points needed to normalize posterior over particle locations
  @inbounds for k in 1:placozoan.observer.nPparticles[]
    try  # sometimes a grid point is out of bounds. Ignored.
      p = placozoan.observer.posterior[gridPoint(placozoan.observer.Pparticle[k])...]
      if p > outlier_threshold
        psum = psum + p
        n = n+1
      end
    catch
    end
  end

  # KLD posterior || particles
  @inbounds for k in 1:placozoan.observer.nPparticles[]
    try  # sometimes a grid point is out of bounds. Ignored.
      p = placozoan.observer.posterior[gridPoint(placozoan.observer.Pparticle[k])...]
      if p > outlier_threshold
        S = S + p*log2(p/psum) # nb divided by psum later to get p/psum*log2(p/psum) 
      end
    catch
    end
  end
  placozoan.observer.KLD[frame] = S/psum + log2(n)

  # KLD of random uniform sample
  S = 0.0
  psum = 0.0
  n = 0
  u = uniformSampleAnnulus(Float64(placozoan.observer.minRange), Float64(placozoan.observer.maxRange), placozoan.observer.nPparticles[])# uniform random sample

  # psum at uniform sample points
  @inbounds for k in 1:placozoan.observer.nPparticles[]
    try
      p = placozoan.observer.posterior[gridPoint(u[k])...]
      if p > outlier_threshold
        psum = psum + p
        n = n + 1
      end
    catch
    end
  end

  # KLD posterior || uniform sample
  @inbounds for k in 1:placozoan.observer.nPparticles[]
    try
      p = placozoan.observer.posterior[gridPoint(u[k])...]
      if p > outlier_threshold
        S = S + p*log2(p/psum) 
      end
    catch
    end
  end
  placozoan.observer.KLD0[frame] = S/psum +  log2(n)


   # KLD of sample from posterior
   S = 0.0
   psum = 0.0
   n = 0
   s = samplePosterior(placozoan)

   @inbounds for k in 1:placozoan.observer.nPparticles[]
    try
      p = placozoan.observer.posterior[gridPoint(s[k])...]
      if p > outlier_threshold
        psum = psum + p
        n = n+1
      end
    catch
    end
  end

   @inbounds for k in 1:placozoan.observer.nPparticles[]
    try
      p = placozoan.observer.posterior[gridPoint(s[k])...]
      if p > outlier_threshold
        S = S + p*log2(p/psum) 
      end
    catch
    end
  end
  placozoan.observer.KLDI[frame] = S/psum +   log2(n)

 # println(  placozoan.observer.KLDI[frame])
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
  Δ = Int64(cround(4*σ))
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

# function sample(D::AbstractArray, N::Int64)
#   # draw sample s of size n from 2D empirical distribution D by rejection]
#   # sum(D)==1.
#   # samples are returned as Int64 nx2 indices of D
#   # samples are returned as Int64 nx2 indices of D
#   # Uses rectangular uniform proposal distribution whose 1/2-width shrinks
#   # to 3x the standard deviation of the posterior as it converges. 

#   s = fill(0, N, 2)   # for samples

#   (peak, ipeak) = findmax(D)  # maximum probability and its location
#   σ = 2^((entropy(D)-4.0942)/2)  # standard deviation of target distribution


#   i = 0              # sample size counter
#   X = axes(D,1)
#   Y = axes(D,2)

#   # sample region is 3 sd each side of peak
#   Δ = Int64(round(3*σ))
#   x0 = max(minimum(X), ipeak[1] - Δ)
#   x1 = min(maximum(X), ipeak[1] + Δ)
#   y0 = max(minimum(Y), ipeak[2] - Δ)
#   y1 = min(maximum(Y), ipeak[2] + Δ) 

#   while i<N
#     x = rand(x0:x1,1)[]    # uniform random point
#     y = rand(y0:y1,1)[]
#     if peak*rand()[] < D[x,y]  # accept/reject
#       i = i + 1
#       s[i,1] = x
#       s[i,2] = y
#     end
#   end

#   s
# end

# random sample from posterior
# sample size = nPparticles
function samplePosterior(placozoan::Placozoan)

  maxPosterior = maximum(placozoan.observer.posterior)
  s = fill(Point2(0.0, 0.0), placozoan.observer.nPparticles[])  
  n = 0
  while n < placozoan.observer.nPparticles[]
    candidate = uniformSampleAnnulus(placozoan.radius, Float64(placozoan.observer.maxRange), 1)  # random point on mat
    if (maxPosterior*rand()) < placozoan.observer.posterior[gridPoint(candidate[])...]   # rejection sampling
      n = n+1
      s[n] = candidate[]   # accept
    end
  end

  return(s)

end

# posterior probability of predator in M-cell RF
function posteriorInMcellRF(p::Placozoan)
   
  # get m-cell vertices
  mc_pts = decompose(Point2f0, Circle(Point2f0(p.mcell.position[][1],p.mcell.position[][2]), p.mcell.r))
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

# uniform random points in annulus 
# i.e. in region bounded by concentric circles with radii (Rinner, Router) 
function uniformSampleAnnulus(Rinner::Float64, Router::Float64, N::Int64)

  r = sqrt.(Rinner^2 .+ rand(N)*(Router^2 - Rinner^2))
  Θ = 2.0*π*rand(N)
  
  Point2.(r.*cos.(Θ), r.*sin.(Θ))

end

# one uniform random point in annulus 
# i.e. in region bounded by concentric circles with radii (Rinner, Router) 
# returns a Point2, not a Vector{Point2} as the other version of this funciton does
function uniformSampleAnnulus(Rinner::Float64, Router::Float64)

  Θ = 2.0*π*rand()
  sqrt(Rinner^2 + rand()*(Router^2 - Rinner^2))*Point2(cos.(Θ), sin.(Θ))

end


# distance from origin to point x
function distance(x::Point2)
  # distance from origin to X
  sqrt(x[1]^2+x[2]^2)
end

# distances from origin to points in vector x
function distance(x::Vector{Point2{Float64}})
  [distance(u) for u in x ]
end

# distance between two points
function distance(x::Point2, y::Point2)
  distance(x-y)
end

# distances between points in vectors p & q
function distance(p::Vector{Point2{Float64}}, q::Vector{Point2{Float64}})
  [distance(p[i]-q[i]) for i in 1:length(p) ]
end

# distance of points in vector q from point p
function distance(p::Point2{Float64}, q::Vector{Point2{Float64}})
  [distance(p-q[i]) for i in 1:length(q) ]
end

# distance of points in vector q from point p
# same as above but with arguments in reverse order so I can be lazy about remembering the order
function distance(q::Vector{Point2{Float64}}, p::Point2{Float64})
  [distance(p-q[i]) for i in 1:length(q) ]
end



# grid point (index into offset array) closest to point
# as tuple (i,j) that can be splatted to array index
function gridPoint(p::Point2{Float64})

  ( Int(cround(p[1])), Int(cround(p[2])) )

end

# rounding function that behaves like c/c++ round
# (n.5 rounds away from zero)
function cround(x::Float64)
  round(x,  RoundNearestTiesAway)
end


# what is the true range when range in map is r?
function projectRangeToWorld(r::Float64, placozoan::Placozoan)

  R = r*(placozoan.observer.maxRange-placozoan.radius)/ 
             placozoan.marginwidth 
end

# what is the range in map when true range is R?
function projectRangeFromWorld(R::Float64, placozoan::Placozoan)

  r = R*placozoan.marginwidth / (placozoan.observer.maxRange-placozoan.radius)

end


function pMap(p::Point2, placozoan::Placozoan)

   d = distance(p)

  r = (d < placozoan.radius + 50.) ? placozoan.observer.collisionRadius[]*sqrt(d/placozoan.radius) : 1.2*placozoan.observer.collisionRadius[]
  
end

function projectMap(p::Point2, placozoan::Placozoan)

  d = distance(p) - placozoan.radius

  #r = ( (distance(p)/placozoan.radius))^.25
  #r = ( (distance(p)-placozoan.radius)/5.0)^.25
  #r = 0.95 + 0.05*(d^2/(160.0 + d^2))
  r = 1.0

end
