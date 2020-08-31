using Makie
#using Unitful

# R = 100.0u"μm"   # body radius
# c =  20.0u"μm"   # cell radius
# D =   0.5u"mm" # max range
# dD =  5.0u"μm"
#
#
# I = 1.5u"nA"  # current per cell
# ρ = 25.0u"Ω*cm"    # resistivity of seawater 25 Ohm.cm
# h = 10u"μm"     # 10 micron dipole separation

R = 500.0   # body radius
c =  20.0   # cell radius
D =   1500.0 # max range
dD =  5.0


I = 1.5  # current per cell
ρ = 25.0    # resistivity of seawater 25 Ohm.cm
h = 10    # 10 micron dipole separation

nD = Int(round(D-R)/dD)
F = fill(0.0, nD)

# sceneD = ustrip(D)
# sceneR = ustrip(R)
# scene = Scene(resolution = (1000, 1000),
#         limits = FRect(-sceneD, -sceneD, 2*sceneD, 2*sceneD))

#body = poly!(scene, decompose(Point2f0, Circle(Point2f0(0,0), sceneR)),
#            color = :lightgray, strokewidth = 1, strokecolor = :black)

  for r in c:c:(R-c)
    n = round(2π*r/c)
    for i in 1:n
      θ = 2π*(i +  rand()/2.)/n
      ρ = c*randn()/4.
      x = (r+ρ)*cos(θ)
      y = (r+ρ)*sin(θ)
      for d in 1:nD
      F[d] = F[d] + ((d+R-x)^2 + y^2)^(-3.0/2.0)

#      scatter!([x], [y], markersize = c/16)
    end
  end
end

# alt display
scene = Scene(); #(limits = FRect(0.0, 0.0, D-R, maximum(F)))
plot!((1:nD).*dD, F)

1+1
display(scene)
