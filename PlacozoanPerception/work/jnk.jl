Δx = 0.1
x = 0:Δx:1
u(x) = sin(2π*x)
using Plots
plot(u,0,1,lw=3)
scatter!(x,u.(x))
plot!(x,u.(x),lw=3)
