module Fig5CurrentDensity

using JLD2
using LinearSolverBenchmarks: Mosfet
using ChargeTransport
using VoronoiFVM
using ExtendableGrids
using LaTeXStrings

# How to use this function?
# 1) includet("post_process/figures/mosfet_2d/fig5_current_density.jl")
# 2) using PythonPlot
# 3) Fig5CurrentDensity.plot_current_density("post_process/figures/mosfet_2d/data_mosfet_2d_ref1.jld2", Plotter=PythonPlot)
function plot_current_density(
        filename;             # "data_mosfet_2d.jld2"
        Plotter = nothing,    # PythonPlot
        save_figure = false,
        title_figure = nothing
    )

    @load filename grid_refinement solution_eq solution biasValues_drain IV_drain


    # Generate corresponding grid
    grid = Mosfet.generate_2Dgrid(grid_refinement)

    #############################################################################
    # Define VoronoiFVM System
    #############################################################################
    parameter_set = Mosfet.tesca_setup
    p = parameter_set()
    constants = ChargeTransport.teSCA_constants

    # Initialize Data instance
    data = Data(grid, p.numberOfCarriers, constants = constants)

    # model
    data.modelType = Stationary

    # statistics
    data.F .= Boltzmann # Tesca Manual S.84 IFERMI = 0

    # recombination
    data.bulkRecombination = set_bulk_recombination(;
        iphin = p.iphin, iphip = p.iphip,
        bulk_recomb_Auger = true,
        bulk_recomb_radiative = false,
        bulk_recomb_SRH = true,
    )

    # boundary model
    data.boundaryType[p.bregion_gate] = GateContact
    data.boundaryType[p.bregion_drain] = OhmicContact
    data.boundaryType[p.bregion_source] = OhmicContact
    data.boundaryType[p.bregion_bulk] = OhmicContact

    # flux discretization - depends on statistic (Boltzmann -> Scharfetter Gummel)
    data.fluxApproximation .= ScharfetterGummel

    # add parameters to data
    data.params = ChargeTransport.Params(p)

    # Definition ChargeTransport System
    ctsys = ChargeTransport.System(grid, data, unknown_storage = :sparse)
    #############################################################################

    nf = -VoronoiFVM.nodeflux(ctsys.fvmsys, solution)
    j = nf[:, p.iphin, :] ./ (1.0e4 * p.zaus) # so that A/ m^2 (zaus in m)

    jx = j[1, :]; jy = j[2, :]
    jabs = sqrt.(jx .^ 2 .+ jy .^ 2)

    # For the arrows
    coord = grid[Coordinates]
    x = coord[1, :]; y = coord[2, :]

    L = 0.03 * (maximum(y) - minimum(y))     # length arrows
    eps = 1.0e-30
    ux = jx ./ (jabs .+ eps);  uy = jy ./ (jabs .+ eps)
    idx = 1:1:length(x)

    Plotter.figure()
    Plotter.tripcolor(
        Mosfet.tridata(grid)...,
        jabs,
        norm = Plotter.matplotlib.colors.LogNorm(vmin = 1.0e-8, vmax = 1.0e4),
        shading = "gouraud",
        rasterized = true
    )
    cbar = Plotter.colorbar(orientation = "vertical", extend = "both")
    Plotter.quiver(
        x[idx], y[idx], (L * ux)[idx], (L * uy)[idx],
        color = "black",
        scale = 1.0,
        angles = "xy",
        scale_units = "xy",
        width = 0.002 # thickness arrows
    )

    Plotter.xlabel(L"\mathrm{width\ [m]}", fontsize = 15)
    Plotter.ylabel(L"\mathrm{height\ [m]}", fontsize = 15)

    ax = Plotter.gca()
    ax.xaxis.set_major_locator(Plotter.matplotlib.ticker.MultipleLocator(0.5e-6))
    ax.yaxis.set_major_locator(Plotter.matplotlib.ticker.MultipleLocator(0.5e-6))
    ax.tick_params(axis = "both", which = "major", labelsize = 10)
    ax.xaxis.get_offset_text().set_fontsize(10)
    ax.yaxis.get_offset_text().set_fontsize(10)
    cbar.ax.tick_params(labelsize = 10)
    cbar.set_label(L"|\mathrm{j}| [\mathrm{A}/\mathrm{cm}^2]", fontsize = 15)

    current_figure = Plotter.gcf()

    if save_figure
        Plotter.savefig(title_figure, bbox_inches = "tight", pad_inches = 0.02)
    end

    return display(current_figure)
end # function

end # module
