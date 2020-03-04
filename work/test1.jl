using Makie
using MakieLayout
using Distributions

scene = Scene(resolution = (600,600), camera=campixel!)
scene_layout = GridLayout(scene, 2, 1,
                    rowsizes = [Relative(0.25), Relative(0.75)],
                    alignmode = Outside(30, 30, 30, 30))

wald_layout = GridLayout(3,1,
         rowsizes = [Relative(.05), Relative(.05), Relative(.9)])

wald_layout[1,1] = tau_slider = LSlider(scene, range=LinRange(-2.0, 2.0, 101))
tau_slider.value[] = 0.0
wald_layout[2,1] = lam_slider = LSlider(scene, range=LinRange(2.0, 4.0, 101))
wald_layout[3,1]= wald_plot_axis = LAxis(scene,
                                        yticksvisible = false,
                                        yticklabelsvisible = false)
exwald_plot_axis.xlabel = "Interval (ms)"
exwald_plot_axis.ylabel = "probability density"

function waldpdf(log_μλ, x)


   μ = 10.0^log_μλ[1]
   λ = 10.0^log_μλ[2]

    return pdf.(InverseGaussian(μ,λ), x)
end


wald_x = collect(0.0:.5:100.0)
# wald_pdf = waldpdf( 12.0,
#                         10.0^lam_slider.value[],
#                         10.0^tau_slider.value[],
#                         exwald_x)

q = Node((tau_slider.value[], lam_slider.value[]))
exwald_pdf_plothandle = plot!(wald_plot_axis, exwald_x,
                  lift(x ->  waldpdf( 12.0, x, exwald_x),  q ) )
                              )
