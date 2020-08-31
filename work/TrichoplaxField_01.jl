using Makie


preyRadius= 400.   # body radius
preyMargin = 100.   # width of margin
c = 20.   # cell radius
D = 5.0e3 # max range (from centre) in μm
dD = 10  # grid spacing for field evaluation

nD = Int(D-preyRadius)

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

# function computes field strength at radius r
# μV/cm
E_(r) = 1.0e6*2π*ρ*I*δ./r.^3
# array holds computed field strength from edge of Trichoplax
E = fill(0.0, nD)

# channel open probability
kT = kB*T
v0 = -σ*log(0.1/(1.0-0.1))
p(e) =  1.0./(1 .+ exp.(-(e.-v0)/σ))
scene = Scene(resolution = (1000, 1000),
limits = FRect(-5*preyRadius, -5*preyRadius, 10*preyRadius, 10*preyRadius))

body = poly!(scene, decompose(Point2f0, Circle(Point2f0(0,0), preyRadius)),
            color = :lightgray, strokewidth = 1, strokecolor = :black)

# fill body with layers of cells
# compute field strength at distance d from edge of body
for a in c:c:(preyRadius-c)
  n = round(2π*a/c)
  x = [ a*cos(2π*i/n) for i in 1:n]
  y = [ a*sin(2π*i/n) for i in 1:n]
  for d in 1:nD
    r = sqrt.(((d.+preyRadius.-x).^2 + y.^2)).*1.0e-4
    E[d] = E[d] + sum(E_(r))
  end
  #scatter!(x, y, markersize = c/16)
end

#display(scene)




# field strength at distance from animal
d = (1:nD).*1.0e-4
plt = plot(d, E)

# voltage μV
V = cumsum(E)*1.0e-4
V = V[end].-V
Vnoisy = V + 1.0e6*σ*randn(size(V))

# plot voltage (μV)
d0 = 1
d1 = 2000
noisyPlot = plot(d[d0:d1]*1.0e4,Vnoisy[d0:d1], color = :lightblue)
voltagePlot =plot!(d[d0:d1]*1.0e4,V[d0:d1], color = :magenta, linewidth = 1)

noisyPlot[Axis][:names][:axisnames] = ("μm from edge", "μV")
openPlot = plot!(d[d0:d1]*1.0e4,.04*p(V[d0:d1]*1e-6), color = :blue)
plot!([d[d0], d[d1]]*1.0e4,.04*p(0)*[1,1], color = :gray)
display(voltagePlot)
title(voltagePlot, "Voltage across receptor")
