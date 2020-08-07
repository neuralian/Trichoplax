using DiffEqOperators
using Makie
using Colors
using ColorSchemes

colormaps = collect(AbstractPlotting.all_gradient_names)

"""
  Diffusion (1 step of Heat Equation)
  ∇u(t) = Δu(x,y)
  on square domain (length(x) = length(y))
"""

function diffuse!(x::Array{Float64,1},
                  y::Array{Float64,1},
                  u::Array{Float64,2},
                  α::Float64,
                  Δt,
                  Δx  )

  v = copy(u)
  @inbounds for i in 2:(length(x)-1), j in 2:(length(y)-1)

    u[i,j]+= Δt*α*(v[i+1,j] - 4.0*v[i,j] + v[i-1,j] + v[i,j+1] + v[i, j-1])/Δx^2

  end

end



Δx = .01
x = collect(Δx:Δx:(1.0-Δx))
y = x
M = length(x)
v = fill(0.0, M, M)

m0 = Int64(round(M/2))
dm = Int64(round(M/10))

for i in (m0-dm):(m0+dm), j in (m0-dm):(m0+dm)
   v[i,j] = 1.0
 end

 D = 100
 limits=FRect(0.0, 0.0, 1.0, 1.0)
 scene = Scene(resolution = (800,800), scale_plot = false,
               show_axis = false, limits=limits)


heatmap!(collect(x),collect(y), v,
               limits = limits,
               colormap = colormaps[64])
handle = scene[end][3]

α = 0.005
Δt = .005

for i in 1:100
 diffuse!(x,y,v,α, Δt, Δx )
 handle[] = v

 # heatmap!(x,y, v,
 #               limits = limits,
 #               colormap = colormaps[64])
 display(scene)
 sleep(.05)

end
