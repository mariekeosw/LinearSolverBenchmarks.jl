module Fig4ComparisonHoleDensity

using JLD2
using ChargeTransport
using LinearSolverBenchmarks: Mosfet
using LaTeXStrings

# How to use this function?
# 1) includet("post_process/figures/mosfet_2d_tesca/fig4_comparison_hole_density.jl")
# 2) using PythonPlot
# 3) Fig4ComparisonHoleDensity.comparison_hole_density("post_process/figures/mosfet_2d_tesca/tesca_files/data_mosfet_ref0_dio.dat", "post_process/figures/mosfet_2d_tesca/ct_files/data_mosfet_2d_ref0.jld2", Plotter=PythonPlot)
function comparison_hole_density(
        datafile_tesca::String, # Solution data Tesca, dio.dat file
        datafile_ct::String;    # Solution data Charge Transport, .jld2 file
        Plotter = nothing,      # PythonPlot
        save_figure = false,
        title_figure = nothing
    )

    # Read tesca data from dio.dat file
    data = ChargeTransport.read_diodat(datafile_tesca)
    np_tesca = data["hDensity"]

    # Read charge transport data from .jld2 file
    @load datafile_ct grid_refinement solution_eq solution biasValues_drain IV_drain

    # Get MOSFET parameter setup and constants
    parameter_set = Mosfet.tesca_setup
    p = parameter_set()
    constants = ChargeTransport.teSCA_constants

    # Calculate hole density
    np_ct = p.Nv .* exp.(p.zp * (constants.q * (solution[p.iphip, :] .- solution[p.ipsi, :]) .+ p.Ev) ./ (constants.k_B * p.T))
    np_ct_cm = np_ct * 1.0e-6

    # Create the underlying grid (used to create the tesca and ct data above,
    # automatically connected to the ct data as it uses the same grid_refinement)
    grid = Mosfet.generate_2Dgrid(grid_refinement)

    # Plot Error
    err = abs.(np_tesca .- np_ct_cm)

    Plotter.figure()
    Plotter.tripcolor(
        Mosfet.tridata(grid)..., err, # change to np_tesca or np_ct_cm here
        norm = Plotter.matplotlib.colors.LogNorm(vmin = 1.0e3, vmax = 1.0e16),
        shading = "gouraud",
        rasterized = true
    )

    cbar = Plotter.colorbar(orientation = "vertical", extend = "both")

    Plotter.xlabel(L"\mathrm{width\ [m]}", fontsize = 15)
    Plotter.ylabel(L"\mathrm{height\ [m]}", fontsize = 15)

    ax = Plotter.gca()
    ax.tick_params(axis = "both", which = "major", labelsize = 10)
    ax.xaxis.get_offset_text().set_fontsize(10)
    ax.yaxis.get_offset_text().set_fontsize(10)
    cbar.ax.tick_params(labelsize = 10)
    cbar.set_label(L"\mathrm{|error|\ [cm^{-3}]}", fontsize = 15)

    current_figure = Plotter.gcf()

    if save_figure
        Plotter.savefig(title_figure, bbox_inches = "tight", pad_inches = 0.02)
    end

    return display(current_figure)

end # function

end # module
