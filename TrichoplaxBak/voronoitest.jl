using VoronoiDelaunay
using Plots
plotly()

using VoronoiDelaunay
tess = DelaunayTessellation()
width = max_coord - min_coord
a = Point2D[Point(min_coord + (i % 10) * width,
            floor(i/10) * width)
            for i in 1:100]
push!(tess, a)

# push!(tess, Point(1.4, 1.4))
# push!(tess, Point(1.6, 1.6))
# push!(tess, Point(1.4, 1.6))
# push!(tess, Point(1.6, 1.4))
# push!(tess, Point(1.5, 1.5))
# push!(tess, Point(1.5, 1.6))
# push!(tess, Point(1.5, 1.4))

xD,yD = getplotxy(delaunayedges(tess))
xV,yV = getplotxy(voronoiedges(tess))

plot(xD,yD, xlims = (1,2), ylims = (1,2),
   aspectratio = 1.0)
plot!(xV,yV)
