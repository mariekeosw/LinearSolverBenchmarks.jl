module Fig6IvCurve

using JLD2
using LaTeXStrings

# How to use this function?
# 1) includet("post_process/figures/mosfet_2d/fig6_iv_curve.jl")
# 2) using PythonPlot
# 3) Fig6IvCurve.plot_iv("post_process/figures/mosfet_2d/data_mosfet_2d_ref1.jld2", Plotter=PythonPlot)
function plot_iv(
        filename;             # "data_mosfet_2d.jld2"
        Plotter = nothing,    # PythonPlot
        save_figure = false,
        title_figure = nothing
    )

    @load filename grid_refinement solution_eq solution biasValues_drain IV_drain

    Plotter.figure()
    Plotter.plot(biasValues_drain, IV_drain, marker = "o", markersize = 3)
    Plotter.grid()

    Plotter.xlabel(L"\mathrm{voltage\ drain\ [V]}", fontsize = 15)
    Plotter.ylabel(L"\mathrm{total\ current\ drain\ [A]}", fontsize = 15)

    ax = Plotter.gca()
    ax.xaxis.set_major_locator(Plotter.matplotlib.ticker.MultipleLocator(1.0))
    ax.tick_params(axis = "both", which = "major", labelsize = 10)

    current_figure = Plotter.gcf()

    if save_figure
        Plotter.savefig(title_figure, bbox_inches = "tight", pad_inches = 0.04)
    end

    return display(current_figure)
end


end # module
