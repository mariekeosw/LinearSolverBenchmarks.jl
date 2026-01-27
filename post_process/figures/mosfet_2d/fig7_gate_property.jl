module Fig7GateProperty

using JLD2
using LaTeXStrings
using ColorSchemes

# How to use this function?
# 1) includet("post_process/figures/mosfet_2d/fig7_gate_property.jl")
# 2) using PythonPlot
# 3) Fig7GateProperty.plot_gate_property("post_process/figures/mosfet_2d/data_mosfet_2d_gate_property.jld2", Plotter=Plots)
function plot_gate_property(filename; Plotter = nothing, save_figure = false, title_figure = nothing)

    if isnothing(Plotter)
        error("Please add Plots.jl in local environment and use Plotter = Plots")
    end

    # load data from .jld2 file
    @load filename maxbias_gate biasValues_drain IV_drain

    cols = ColorSchemes.seaborn_deep

    p = Plotter.plot(
        xlabel = L"\mathrm{voltage\ drain\ [V]}",
        ylabel = L"\mathrm{total\ current\ drain\ [A]}",
        title = "",
        legend = :topleft,
        bottom_margin = 10 * Plotter.px,
        guidefontsize = 16,
        framestyle = :box,
        foreground_color_border = :black,
        linewidth = 2,
        grid = true,
        gridlinewidth = 0.15,
        gridcolor = :black,
        gridalpha = 0.6,
    )

    for (i, gate_bias) in enumerate(maxbias_gate)
        Plotter.plot!(
            p,
            biasValues_drain[gate_bias],
            IV_drain[gate_bias];
            linecolor = cols[i],
            label = "voltage gate: $(gate_bias) V",
            linestyle = :solid,
            marker = :circle,
            markersize = 3,
            markerstrokewidth = 0,
            markercolor = cols[i],
        )
    end

    if save_figure
        Plotter.savefig(p, title_figure)
    end

    return p
end # plot_gate_property

end # module
