# Cell structure
# dev for Trichoplax model
# Mike Paulin University of Otago 2019
#


using Makie
using AbstractPlotting
using Colors
using Distributions

#
cellCount = Int32(0)


# Cell structure
struct Cell2D
  who:: Int32
  x:: Array{Float64,1}
  y:: Array{Float64,1}
  r   # edge rest length
  k   # edge stiffness
  σ   # edge energy
end

# Cell constructor
# nb faces indexed anticlockwise from 1 ↗
Cell2D(x₀, y₀, r) = Cell2D( cellCount,
                          x₀ .+ [r*cos(i*π/3) for i in 0:6],
                          y₀ .+ [r*sin(i*π/3) for i in 0:6],
                          r, 1., 1.)


# # Cell sprouter
# # sprout a cell from existing cell
function sprout(cell, whatface)

  x0 = cell.x[1]
  y0 = cell.y[1]
  dx = cell.r*cos(π/3)
  dy = cell.r*sin(π/3)

  #Cell2D(cell.x[2]  + 1.0, cell.y[2], cell.r)
  #Cell2D(cell.x[1] - 1.0, cell.y[1] + 2.0*dy, cell.r)
  #Cell2D(cell.x[3] - 1.0, cell.y[3], cell.r)
  #Cell2D(cell.x[5] -1.0, cell.y[5], cell.r)
  #Cell2D(cell.x[1] -1.0, cell.y[1] - 2.0*dy, cell.r)
  Cell2D(cell.x[6] + 1.0, cell.y[6], cell.r)
end

# Cell renderer
function plot(cell)
  # draw cell
  poly!(Point2f0[[cell.x[i], cell.y[i]] for i in 1:6],
     color = RGBA(1.0, 0.5, 0.5, .75), strokecolor = :black,
     strokewidth = 1, scale_plot = false,
     backgroundcolor = RGBA(0.0, 0.5, 0.0, 0.25) )
end

Mat = Scene(limits = FRect(-10, -10, 20, 20))

cellCount +=1
cell = Cell2D(0,0,1)
plot(cell)

cellCount +=1
cell = sprout(cell, 1)
plot(cell)

# draw  cell

print("hello")

display(Mat)

#sprout(1,1)
