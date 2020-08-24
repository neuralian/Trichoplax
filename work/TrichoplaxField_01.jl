using Makie


preyRadius = 200.   # body radius
preyMargin = 100.   # width of margin
c = 20.   # cell radius
D = 1.0e4 # max range (from centre)
dD = 10.  # grid spacing for field evaluation

nD = Int(D-preyRadius)
V = fill(0.0, nD)

scene = Scene(resolution = (1000, 1000),
limits = FRect(-5*preyRadius, -5*preyRadius, 10*preyRadius, 10*preyRadius))

body = poly!(scene, decompose(Point2f0, Circle(Point2f0(0,0), preyRadius)),
            color = :lightgray, strokewidth = 1, strokecolor = :black)

# fill body with layers of cells
# compute field strength at distance d from edge of body
for r in c:c:(preyRadius-c)
  n = round(2π*r/c)
  x = [ r*cos(2π*i/n) for i in 1:n]
  y = [ r*sin(2π*i/n) for i in 1:n]
  for d in 1:nD
    V[d] = V[d] + sum(((d.+preyRadius.-x).^2 + y.^2).^-1.5)
  end
  scatter!(x, y, markersize = c/16)
end
#plot(1:nD, F)
display(scene)

# Johnson-Nyquist noise
kB = 1.38e-23
T = 300.
Ω = 20.e6   # 20MΩ
Δf = 1.0e3   # 1KHz
σ = sqrt(4.0*kB*T*Δf)
print(σ)
