using Makie


R = 500.   # body radius
c = 20.   # cell radius
D = 1000. # max range
dD = 1.

nD = Int(D-R)
F = fill(0.0, nD)

scene = Scene(resolution = (1000, 1000), limits = FRect(-2*R, -2*R, 4*R, 4*R))

body = poly!(scene, decompose(Point2f0, Circle(Point2f0(0,0), R)),
            color = :lightgray, strokewidth = 1, strokecolor = :black)

  for r in c:c:(R-c)
    n = round(2π*r/c)
    for i in 1:n
      θ = 2π*(i +  rand()/2.)/n
      ρ = c*randn()/4.
      x = (r+ρ)*cos(θ)
      y = (r+ρ)*sin(θ)
      for d in 1:nD
      F[d] = F[d] + ((d+R-x)^2 + y^2)^(-3.0/2.0)

      #scatter!([], [], markersize = c/16)
    end
  end
end
plot(1:nD, F)
#display(scene)
