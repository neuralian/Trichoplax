# Placozoan_Dev.jl

# function trichoplax_sim()

# MAIN

bodylayers = 8 # number of body cell layers
margin = 3  # number of layers in gut margin ("brain")
range = 1.25   # sensory range in body radius units
celldiameter = 10.0
skeleton_springconstant= 5.0e-2
cell_pressureconstant = 1.0e0
cell_surface_energy_density  = 5.0e1
dt = .001

param = trichoplaxparameters(   bodylayers,
                                margin,
                                range,
                                skeleton_springconstant,
                                cell_pressureconstant,
                                cell_surface_energy_density,
                                celldiameter,
                                dt)

@time trichoplax = Trichoplax(param)
# trichoplax.param.k2[] = 5.0e-2    # cytoskeleton spring constant /2
# trichoplax.param.σ[]  = 5.0e1   # cell surface energy density
# trichoplax.param.ρ[]  = 1.0e0 #1.0e2    # cell turgor pressure energy per unit volume


# Draw
R = bodylayers*celldiameter    # approx radius of Trichoplax (for scene setting)
D = Int64(3*R)  # scene diameter
limits=FRect(-D, -D, 2*D, 2*D)
scene = Scene(resolution = (800,800), scale_plot = false,
              show_axis = false, limits=limits)


# electric field
efieldcolormap = ColorScheme(
                [RGB{Float64}(
                    (i<.85) ? .25 : i,
                    0.8-i,
                    i
                    ) for i in 0:0.1:1.0])

x = y = collect(-D:4:D)
m = length(x)
n = length(y)
algamat = fill(0.0, m,n)

for i in 1:m
    for j in 1:n
        algamat[i,j] = sqrt(x[i]^2 + y[j]^2)/sqrt(D^2 + D^2)
    end
end

heatmap!(x,y, algamat,
                limits = limits,
                colormap = efieldcolormap)

# # scatter bacteria (point objects) over the scene
# nbacteria = 50
# bactrect = FRect(20,20, 30, 30)
# bacteria = growbacteria(nbacteria, bactrect)

# draw trichoplax cells
cells_handle = draw(scene, trichoplax, RGB(.25, .25, .25), 1)

# colour the cells

imagecells( trichoplax,
            1:trichoplax.anatomy.nstomach,
            trichoplax.state.potential, 1)
imagecells( trichoplax,
            (trichoplax.anatomy.nstomach+1):trichoplax.anatomy.ncells,
            trichoplax.state.potential, 2)
# ch = potentialmap(scene, trichoplax)

#
reflect(trichoplax)

display(scene)

# project
# bodyCenter = colmeans(verticesofcell(1, trichoplax))
# skinvertex = getskinvertexcoords(trichoplax)
# bodyRadius = meanvec(distance(bodyCenter, skinvertex)[:])
# plotcircle!(scene, bodyCenter, bodyRadius, color = :green)
# viewRadius = 1.5*bodyRadius
# plotcircle!(scene, bodyCenter, viewRadius, color = :red)
# gutboundaryvertex = getgutboundaryvertexcoords(trichoplax)
# gutRadius = meanvec(distance(bodyCenter, gutboundaryvertex)[:])
# plotcircle!(scene, bodyCenter, gutRadius, color = :green)





# restvolume = copy(trichoplax.state.volume)
# i0 = 4
# i1 = vcat(i0, trichoplax.anatomy.neighbourcell[i0,:])

# # record(scene, "trichoplaxdev.mp4", 1:100) do tick
# for tick in 1:100
#     global trichoplax
#     global scene
#     if tick < 50
#         whichbacteriahere = bacteriahere(bacteria, trichoplax)
#         for j in 1:trichoplax.anatomy.stomach
#             for k in 1:length(whichbacteriahere[j])
#                 trichoplax.state.potential[j] = trichoplax.state.potential[j] -
#                         1.0/(bacteria.deadticks[whichbacteriahere[j][k]] .+ 1)
#                 bacteria.deadticks[whichbacteriahere[j][k]] = bacteria.deadticks[whichbacteriahere[j][k]] .+ 1
#             end
#         end
#     end
#
#
#     trichoplax = diffusepotential(trichoplax,600)
#
#     trichoplax.state.volume[:] =
#         restvolume.*(1.0 .-
#         sign.(trichoplax.state.potential).*
#         sqrt.(abs.(trichoplax.state.potential/4.0)))
#     trichoplax = morph(trichoplax, .0001, 25)
#
#     redraw(trichoplax,cells_handle)
#     potential_remap(trichoplax, ch, 1)
#
#     println(tick)
#     display(scene)
#     sleep(.005)
# end

# end


# ch[12][:color] = [   RGB{Float64}(0.913603,0.0,0.0),
#                                 RGB{Float64}(0.,0.,1.0),
#                                 RGB{Float64}(0.913603,0.0,0.788739),
#                                 RGB{Float64}(0.913603,0.61886,0.788739),
#                                 RGB{Float64}(0.913603,0.61886,0.788739),
#                                 RGB{Float64}(0.913603,0.61886,0.788739),
#                                 RGB{Float64}(0.913603,0.61886,0.788739)]
