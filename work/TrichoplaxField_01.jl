using Makie


preyRadius = 400.   # body radius
preyMargin = 100.   # width of margin
c = 20.   # cell radius
D = 5.0e3 # max range (from centre) in μm
dD = 10  # grid spacing for field evaluation

nD = Int(D-preyRadius)

# source
ρ = 25.0   # Resisitivity of seawater 25Ω.cm
δ = 20.e-6*100.  # dipole separation 10μm in cm
I = 5.e-10     # dipole current 0.5nA

# function computes field strength at radius r
# μV/cm
E_(r) = 1.0e6*2π*ρ*I*δ./r.^3
# array holds computed field strength from edge of Trichoplax
E = fill(0.0, nD)

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
    r = sqrt.(sum(((d.+preyRadius.-x).^2 + y.^2))).*1.0e-4
    E[d] = E[d] + E_(r)
  end
  #scatter!(x, y, markersize = c/16)
end

#display(scene)

# Johnson-Nyquist noise
kB = 1.38e-23
T = 300.
Ω = 20.e6   # 20MΩ
Δf = 1.0e3   # 1KHz
σ = sqrt(4.0*kB*T*Δf)
print(σ)


# field strength at distance from animal
d = (1:nD).*1.0e-4
plt = plot(d, E*1.e4)

# voltage μV
V = cumsum(E)
V = V[end].-V
# plot(d,V)

display(plt)
