module Fig4HoleDensity

using JLD2
using LinearSolverBenchmarks: Mosfet
using ChargeTransport
using ExtendableGrids
using LaTeXStrings

# How to use this function?
# 1) includet("post_process/figures/mosfet_2d/fig4_hole_density.jl")
# 2) using PythonPlot
# 3) Fig4HoleDensity.plot_hole_density("post_process/figures/mosfet_2d/data_mosfet_2d_ref1.jld2", Plotter=PythonPlot)
function plot_hole_density(
        filename;             # "data_mosfet_2d.jld2"
        Plotter = nothing,    # PythonPlot
        Equilibrium = false,
        save_figure = false,
        title_figure = nothing
    )

    @load filename grid_refinement solution_eq solution biasValues_drain IV_drain

    # Choose solution in equilibrium if Equilibrium = true
    if Equilibrium
        solution = solution_eq
    end

    # Get MOSFET parameter setup and constants
    parameter_set = Mosfet.tesca_setup
    p = parameter_set()
    constants = ChargeTransport.teSCA_constants

    # Calculate hole density
    np = p.Nv .* exp.(p.zp * (constants.q * (solution[p.iphip, :] .- solution[p.ipsi, :]) .+ p.Ev) ./ (constants.k_B * p.T))
    np_cm = np * 1.0e-6

    # Generate corresponding grid
    grid = Mosfet.generate_2Dgrid(grid_refinement)

    Plotter.figure()
    Plotter.tripcolor(
        Mosfet.tridata(grid)..., np_cm,
        norm = Plotter.matplotlib.colors.LogNorm(vmin = 1.0e3, vmax = 1.0e16),
        shading = "gouraud",
        rasterized = true
    )
    cbar = Plotter.colorbar(orientation = "vertical", extend = "both")

    Plotter.xlabel(L"\mathrm{width\ [m]}", fontsize = 15)
    Plotter.ylabel(L"\mathrm{height\ [m]}", fontsize = 15)

    ax = Plotter.gca()
    ax.xaxis.set_major_locator(Plotter.matplotlib.ticker.MultipleLocator(0.5e-6))
    ax.yaxis.set_major_locator(Plotter.matplotlib.ticker.MultipleLocator(0.5e-6))
    ax.tick_params(axis = "both", which = "major", labelsize = 10)
    ax.xaxis.get_offset_text().set_fontsize(10)
    ax.yaxis.get_offset_text().set_fontsize(10)
    cbar.ax.tick_params(labelsize = 10)
    cbar.set_label(L"\mathrm{density}\ [\mathrm{cm}^{-3}]", fontsize = 15)

    current_figure = Plotter.gcf()

    if save_figure
        Plotter.savefig(title_figure, bbox_inches = "tight", pad_inches = 0.02)
    end

    return display(current_figure)

end # function

end # module
