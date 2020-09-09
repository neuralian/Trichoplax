using Plots
pyplot()

kT = 1.38e-23 * 293.
λ = kT/10.
α = 2.5
ro = 1.0
R = 0.5
Fo = λ * ro ^ -α
r = 0.001:0.001:2.0

Ez(r) = λ*(r .+ R).^-α


Po(r) = 1.0 ./(1.0 .+ exp.((Ez(ro) .- Ez(r))/kT))

plot(r, Po(r))
plot!(r, Ez(r))
