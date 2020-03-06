# Placozoan_Dev.jl

# function trichoplax_sim()

# MAIN
<<<<<<< Updated upstream
=======
bodylayers = 4 # number of body cell layers
# mapdepth = 1     # map layers
celldiam = 10.0
SceneWidth = bodylayers*celldiam
>>>>>>> Stashed changes


bodylayers = 3 # number of body cell layers
celldiameter = 10.0
skeleton_springconstant= 5.0e-2
cell_pressureconstant = 1.0e0
cell_surface_energy_density  = 5.0e1
dt = .001

param = trichoplaxparameters(   bodylayers,
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
D = 3*R  # scene diameter
limits=FRect(-D/2, -D/2, D, D)
scene = Scene(resolution = (800,800), scale_plot = false,
              show_axis = false, limits=limits)

# scatter bacteria (point objects) over the scene
nbacteria = 100
bacteria = growbacteria(nbacteria, limits)

# draw trichoplax cells
cells_handle = draw(scene, trichoplax, RGB(.25, .25, .25), 1)

# colour the cells
ch = potentialmap(scene, trichoplax)
<<<<<<< Updated upstream

display(scene)

restvolume = copy(trichoplax.state.volume)
i0 = 4
i1 = vcat(i0, trichoplax.anatomy.neighbourcell[i0,:])

#record(scene, "trichoplaxdev.mp4", 1:200) do i
for i in 1:25
=======
# display(scene)
# sleep(0.25)
#
restvolume = copy(trichoplax.volume)
i0 = 8
i1 = vcat(i0, trichoplax.neighbourcell[i0,:])

#record(scene, "trichoplaxdev.mp4", 1:200) do i
for i in 1:5
>>>>>>> Stashed changes
    global trichoplax
    global scene
    q1 = i<=50 ? 1 : 0
    q2 = i<=100 ? 1 : 0
    trichoplax.state.potential[i1] .= 2.0*q1 - 0.025*(1-q1)*q2
    # trichoplax.potential[9] = i<=50 ? 1.5 : 0.0

    trichoplax = diffusepotential(trichoplax,600)

    trichoplax.state.volume[:] =
        restvolume.*(1.0 .-
        sign.(trichoplax.state.potential).*
        sqrt.(abs.(trichoplax.state.potential/4.0)))
    trichoplax = morph(trichoplax, .0001, 25)
    # scene = Scene(resolution = (800,800), scale_plot = false,
    #     limits=FRect(-SceneWidth, -SceneWidth, 2*SceneWidth, 2*SceneWidth))
    redraw(trichoplax,cells_handle)
    potential_remap(trichoplax, ch, 1)
    # potentialmap(scene, trichoplax)
    display(scene)
    sleep(.005)
end

# end


# ch[12][:color] = [   RGB{Float64}(0.913603,0.0,0.0),
#                                 RGB{Float64}(0.,0.,1.0),
#                                 RGB{Float64}(0.913603,0.0,0.788739),
#                                 RGB{Float64}(0.913603,0.61886,0.788739),
#                                 RGB{Float64}(0.913603,0.61886,0.788739),
#                                 RGB{Float64}(0.913603,0.61886,0.788739),
#                                 RGB{Float64}(0.913603,0.61886,0.788739)]
