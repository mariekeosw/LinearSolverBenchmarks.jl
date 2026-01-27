module Fig2ComparisonElectrostaticPotential

using JLD2
using ChargeTransport
using LinearSolverBenchmarks: Mosfet
using LaTeXStrings

# How to use this function?
# 1) includet("post_process/figures/mosfet_2d_tesca/fig2_comparison_electrostatic_potential.jl")
# 2) using PythonPlot
# 3) Fig2ComparisonElectrostaticPotential.comparison_electrostatic_potential("post_process/figures/mosfet_2d_tesca/tesca_files/data_mosfet_ref0_dio.dat", "post_process/figures/mosfet_2d_tesca/ct_files/data_mosfet_2d_ref0.jld2", Plotter=PythonPlot)
function comparison_electrostatic_potential(
        datafile_tesca::String, # Solution data Tesca, dio.dat file
        datafile_ct::String;    # Solution data Charge Transport, .jld2 file
        Plotter = nothing,      # PythonPlot
        save_figure = false,
        title_figure = nothing
    )

    # Read tesca data from diodat file
    data = ChargeTransport.read_diodat(datafile_tesca)
    solution_psi_tesca = data["ElectrostaticPotential"]

    # Read charge transport data from .
    @load datafile_ct grid_refinement solution_eq solution biasValues_drain IV_drain
    ipsi = 3
    solution_psi_ct = solution[ipsi, :]

    # Create the underlying grid (used to create the tesca and ct data above,
    # automatically connected to the ct data as it uses the same grid_refinement)
    grid = Mosfet.generate_2Dgrid(grid_refinement)

    # Plot Error
    err = abs.(solution_psi_tesca .- solution_psi_ct)

    Plotter.figure()
    Plotter.tripcolor(Mosfet.tridata(grid)..., err) # change to solution_psi_tesca or solution_psi_ct here
    cbar = Plotter.colorbar(orientation = "vertical", extend = "both")

    Plotter.xlabel(L"\mathrm{width\ [m]}", fontsize = 15)
    Plotter.ylabel(L"\mathrm{height\ [m]}", fontsize = 15)

    ax = Plotter.gca()
    ax.tick_params(axis = "both", which = "major", labelsize = 10)
    ax.xaxis.get_offset_text().set_fontsize(10)
    ax.yaxis.get_offset_text().set_fontsize(10)
    cbar.ax.tick_params(labelsize = 10)
    cbar.set_label(L"\mathrm{|error|\ [V]}", fontsize = 15)

    #Plotter.tight_layout()
    current_figure = Plotter.gcf()

    if save_figure
        Plotter.savefig(title_figure, bbox_inches = "tight", pad_inches = 0.02)
    end

    return display(current_figure)

end # function

end # module
