using Placozoan
using Makie

# Placozoan package must be on load path.
# In Atom:
#    right-click on folder, select "work in folder"
#    push!(LOAD_PATH, pwd())

diameter = 32
trichoplax = make_trichoplax(diameter);


# Set scene
sceneWidth = Int64(round(diameter*1.5))
if isodd(sceneWidth)
    sceneWidth = sceneWidth + 1
end

world = Scene(limits = FRect(-sceneWidth/2, -sceneWidth/2,
                              sceneWidth,sceneWidth), scale_plot = false)


draw(trichoplax)
display(world)
