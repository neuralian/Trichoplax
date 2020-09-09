using Makie

x = range(0, stop = 2pi, length = 80)
y1 = sin.(x)
y2 = exp.(-x) .* cos.(2pi * x)

scene = lines(x, y1, color = :blue)
scatter!(scene, x, y1, color = :red, markersize = 0.1)

lines!(scene, x, y2, color = :black)
scatter!(scene, x, y2, color = :green, marker = :utriangle, markersize = 0.1)
