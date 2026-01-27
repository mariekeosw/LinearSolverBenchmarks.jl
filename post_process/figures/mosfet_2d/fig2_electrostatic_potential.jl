module Fig2ElectrostaticPotential

using JLD2
using LinearSolverBenchmarks: Mosfet
using ExtendableGrids
using LaTeXStrings

# How to use this function?
# 1) includet("post_process/figures/mosfet_2d/fig2_electrostatic_potential.jl")
# 2) using PythonPlot
# 3) Fig2ElectrostaticPotential.plot_psi("post_process/figures/mosfet_2d/data_mosfet_2d_ref1.jld2", Plotter=PythonPlot)
function plot_psi(
        filename;             # "data_mosfet_2d.jld2"
        Plotter = nothing,    # PythonPlot
        Equilibrium = false,
        save_figure = false,
        title_figure = nothing
    )

    @load filename grid_refinement solution_eq solution biasValues_drain IV_drain

    # Choose solution in equilibrium when Equilibrium = true
    if Equilibrium
        solution = solution_eq
    end

    # Get MOSFET parameter setup
    parameter_set = Mosfet.tesca_setup
    p = parameter_set()

    # Generate corresponding grid
    grid = Mosfet.generate_2Dgrid(grid_refinement)
    XX = grid[Coordinates][1, :]; YY = grid[Coordinates][2, :]

    fig = Plotter.figure()
    ax = fig.add_subplot(111, projection = "3d")
    ax.plot_trisurf(XX[:], YY[:], solution[p.ipsi, :])

    Plotter.xlabel(L"\mathrm{width\ [m]}", fontsize = 13)
    Plotter.ylabel(L"\mathrm{height\ [m]}", fontsize = 13)
    Plotter.zlabel(L"\mathrm{potential\ [V]}", fontsize = 13)
    ax = Plotter.gca()
    ax.tick_params(axis = "x", labelsize = 10)
    ax.tick_params(axis = "y", labelsize = 10)
    ax.tick_params(axis = "z", labelsize = 10)
    ax.xaxis.get_offset_text().set_fontsize(10)
    ax.yaxis.get_offset_text().set_fontsize(10)
    ax.zaxis.get_offset_text().set_fontsize(10)

    #Plotter.tight_layout()
    current_figure = Plotter.gcf()

    if save_figure
        Plotter.savefig(title_figure, bbox_inches = "tight", bbox_extra_artists = [ax.zaxis.label], pad_inches = 0.2)
    end

    return display(current_figure)
end # function

end # module
