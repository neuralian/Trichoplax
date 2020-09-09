# Placozoan_Dev.jl



# MAIN
bodylayers = 4 # number of body cell layers
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
draw(trichoplax, RGB(.25, .25, .25), 1)

# map potential to cell colour
# default colormap = mint; other options in function potentialmap()
potentialmap(trichoplax)
display(scene)
sleep(0.25)

restvolume = copy(trichoplax.volume)
record(scene, "trichoplaxdev.mp4", 1:150) do i
    global trichoplax
    global scene
    display(scene)
    trichoplax.potential[8] = i<=50 ? 2.0 : 0.0
    trichoplax = diffusepotential(trichoplax,100)
    trichoplax.volume[:] = restvolume.*(1.0 .- sqrt.(trichoplax.potential/4.0))
    morph(trichoplax, .0001, 25)
    scene = Scene(resolution = (800,800), scale_plot = false,
        limits=FRect(-SceneWidth, -SceneWidth, 2*SceneWidth, 2*SceneWidth))
    draw(trichoplax,RGB(.25, .25, .25), 1)
    potentialmap(trichoplax)
    sleep(0.1)
    println(i)
end
