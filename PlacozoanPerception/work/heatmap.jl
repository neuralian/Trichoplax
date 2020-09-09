using Makie

colormaps = collect(AbstractPlotting.all_gradient_names)

scene = Scene(resolution = (800,800))

i = i+1
heatmap!(rand(16,16), colormap = colormaps[i])
println(i)
