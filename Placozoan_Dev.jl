# Placozoan_Dev.jl

# function trichoplax_sim()

# MAIN
bodylayers = 8 # number of body cell layers
# mapdepth = 1     # map layers
celldiam = 10.0
SceneWidth = bodylayers*celldiam



@time trichoplax = Trichoplax(bodylayers, celldiam)
trichoplax.k2[] = 5.0e-2    # cytoskeleton spring constant /2
trichoplax.σ[]  = 5.0e1   # cell surface energy density
trichoplax.ρ[]  = 1.0e0 #1.0e2    # cell turgor pressure energy per unit volume


# Draw
scene = Scene(resolution = (800,800), scale_plot = false,
    limits=FRect(-SceneWidth, -SceneWidth, 2*SceneWidth, 2*SceneWidth))
cells_handle = draw(scene, trichoplax, RGB(.25, .25, .25), 1)

# # map potential to cell colour
# # default colormap = mint; other options in function potentialmap()
ch = potentialmap(scene, trichoplax)
# display(scene)
# sleep(0.25)
#
restvolume = copy(trichoplax.volume)
i0 = 48
i1 = vcat(i0, trichoplax.neighbourcell[i0,:])

record(scene, "trichoplaxdev.mp4", 1:200) do i
#for i in 1:150
    global trichoplax
    global scene
    q1 = i<=50 ? 1 : 0
    q2 = i<=100 ? 1 : 0
    trichoplax.potential[i1] .= 2.0*q1 - 0.025*(1-q1)*q2
    # trichoplax.potential[9] = i<=50 ? 1.5 : 0.0

    trichoplax = diffusepotential(trichoplax,600)

    trichoplax.volume[:] = restvolume.*(1.0 .- sign.(trichoplax.potential).*sqrt.(abs.(trichoplax.potential/4.0)))
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
