using Makie


sourceRadius= 400.   # body radius
sourceMargin = 100.   # width of margin
cellDiam = 20.   # cell radius
maxRange = 5.0e3 # max range (from centre) in μm
gridSize = 10  # grid spacing for field evaluation

nD = Int(maxRange-sourceRadius)

# source
ρ = 25.0   # Resisitivity of seawater 25Ω.cm
δ = 20.e-6*100.  # dipole separation 10μm in cm
I = 0.1e-12     # dipole current 0.1pA

# Johnson-Nyquist noise
kB = 1.38e-23
T = 300.
Ω = 20.e6   # 20MΩ
Δf = 1.0e3   # 1KHz
σ = sqrt(4.0*kB*T*Δf)
print(σ)

# function computes field strength at radius r from dipole source
# μV/cm
E_(r) = 1.0e6*2π*ρ*I*δ./r.^3
# array holds computed field strength from edge of Trichoplax
E = fill(0.0, nD)

# channel open probability
kT = kB*T
v0 = -σ*log(0.1/(1.0-0.1))
p(e) =  1.0./(1 .+ exp.(-(e.-v0)/σ))
scene = Scene(resolution = (1000, 1000),
limits = FRect(-5*sourceRadius, -5*sourceRadius, 10*sourceRadius, 10*sourceRadius))

body = poly!(scene, decompose(Point2f0, Circle(Point2f0(0,0), sourceRadius)),
            color = :lightgray, strokewidth = 1, strokecolor = :black)

# fill body with layers of cells
# compute field strength at distance d from edge of body
function fieldStrength(sourceRadius, sourceMargin, nD, E)
  for a in cellDiam:cellDiam:(sourceRadius-cellDiam)
    n = round(2π*a/cellDiam)
    x = [ a*cos(2π*i/n) for i in 1:n]
    y = [ a*sin(2π*i/n) for i in 1:n]
    for d in 1:nD
      r = sqrt.(((d.+sourceRadius.-x).^2 + y.^2)).*1.0e-4
      E[d] = E[d] + sum(E_(r))
    end
    #scatter!(x, y, markersize = cellDiam/16)
  end
  return E
end

E = fieldStrength(sourceRadius, sourceMargin, nD, E)
#display(scene)

function pOpen(d, V)
   i = Int(round(d))
   return i < length(V)-1 ? p(V[i+1]*1.0e-6) : 0.0
 end


# field strength at distance from animal
d = (1:nD).*1.0e-4
#plt = lines(d, E, strokewidth = 0.1, markersize = 0)

# voltage μV
V = cumsum(E)*1.0e-4
V = V[end].-V
Vnoisy = V + 1.0e6*σ*randn(size(V))

# plot voltage (μV)
d0 = 1
d1 = 2000
noisyPlot = lines(d[d0:d1]*1.0e4,Vnoisy[d0:d1], color = :lightblue, linewidth=1)
voltagePlot =lines!(d[d0:d1]*1.0e4,V[d0:d1], color = :magenta, linewidth = 1)

noisyPlot[Axis][:names][:axisnames] = ("μm from edge", "μV")
openPlot = lines!(d[d0:d1]*1.0e4,.04*p(V[d0:d1]*1e-6), color = :blue, linewidth=1)
lines!([d[d0], d[d1]]*1.0e4,.04*p(0)*[1,1], color = :red, linewidth=1)
display(voltagePlot)
title(voltagePlot, "Voltage across receptor")
