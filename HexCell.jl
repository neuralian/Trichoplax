# Draw hexagonal cell
# dev for Trichoplax model
# Mike Paulin University of Otago 2019
#


using Makie
using AbstractPlotting
using Colors
using Distributions



Mat = Scene()

function drawHexCell(mat, x0,y0, L)
  # x0,y0 = center location, L = side length

  c = cos(π/3.0)
  s = sin(π/3.0)
  lines!(x0 .+ [1.,  c, -c, -1.0, -c, c, 1.0],
         y0 .+ [0.0, s,  s,  0.0, -s,  -s, 0.0],
         limits = FRect(-1, -1, 2, 2), scale_plot=false)

  mat[end]  # return handle to cell

end

# draw  cell
Cell_handle = drawHexCell(Mat, 0.0, 0.0, 1)


Mat
