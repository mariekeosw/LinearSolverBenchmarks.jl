module Fig3ElectronDensity

using JLD2
using LinearSolverBenchmarks: Mosfet
using ChargeTransport
using ExtendableGrids
using LaTeXStrings

# How to use this function?
# 1) includet("post_process/figures/mosfet_2d/fig3_electron_density.jl")
# 2) using PythonPlot
# 3) Fig3ElectronDensity.plot_electron_density("post_process/figures/mosfet_2d/data_mosfet_2d_ref1.jld2", Plotter=PythonPlot)
function plot_electron_density(
        filename;             # "data_mosfet_2d.jld2"
        Plotter = nothing,    # PythonPlot
        Equilibrium = false,
        save_figure = false,
        title_figure = nothing
    )

    @load filename grid_refinement solution_eq solution biasValues_drain IV_drain

    # Choose solution in equilibrium wenn Equilibrium = true
    if Equilibrium
        solution = solution_eq
    end

    # Get MOSFET parameter setup and constants
    parameter_set = Mosfet.tesca_setup
    p = parameter_set()
    constants = ChargeTransport.teSCA_constants

    # Calculate electron density
    nn = p.Nc .* exp.(p.zn * (constants.q * (solution[p.iphin, :] .- solution[p.ipsi, :]) .+ p.Ec) ./ (constants.k_B * p.T))
    nn_cm = nn * 1.0e-6

    # Generate corresponding grid
    grid = Mosfet.generate_2Dgrid(grid_refinement)

    Plotter.figure()
    Plotter.tripcolor(
        Mosfet.tridata(grid)..., nn_cm,
        norm = Plotter.matplotlib.colors.LogNorm(vmin = 1.0e6, vmax = 1.0e20),
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
