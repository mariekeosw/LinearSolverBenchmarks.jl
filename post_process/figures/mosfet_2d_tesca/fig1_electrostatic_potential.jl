module Fig1ElectrostaticPotential

using ChargeTransport
using ExtendableGrids
using LaTeXStrings

# Gridfile needs to be in .dom format, dom file in \mu is transformed to cm

# How to use this function?
# 1) includet("post_process/figures/mosfet_2d_tesca/fig1_electrostatic_potential.jl")
# 2) using PythonPlot
# 3) Fig1ElectrostaticPotential.plot_tesca_electrostaticpotential("post_process/figures/mosfet_2d_tesca/tesca_files/grid_ref0.dom", "post_process/figures/mosfet_2d_tesca/tesca_files/data_mosfet_ref0_dio.dat", Plotter = PythonPlot)
function plot_tesca_electrostaticpotential(
        gridfile::String,
        datafile::String;
        Plotter = nothing,
        save_figure = false,
        title_figure = nothing
    )

    Plotter.pygui(false) # Plots as Pop-Ups
    Plotter.close("all")

    # Read grid from dom file
    grid = simplexgrid(gridfile, format = "dom")
    grid[Coordinates] *= 1.0e-6
    trim!(grid)

    # Grid for plotting
    XX = grid[Coordinates][1, :]; YY = grid[Coordinates][2, :]

    # Read data from diodat file
    data = ChargeTransport.read_diodat(datafile)

    # Read Electrostatic Potential from data
    solution_psi = data["ElectrostaticPotential"]

    # For the plotting routines
    ipsi = 1 # solution contains only psi
    solution_psi = data["ElectrostaticPotential"]

    fig = Plotter.figure()
    ax = fig.add_subplot(111, projection = "3d")
    ax.plot_trisurf(XX[:], YY[:], solution_psi)

    Plotter.xlabel(L"\mathrm{width\ [m]}", fontsize = 12)
    Plotter.ylabel(L"\mathrm{height\ [m]}", fontsize = 12)
    Plotter.zlabel(L"\mathrm{potential\ [V]}", fontsize = 12)
    ax = Plotter.gca()
    ax.tick_params(axis = "x", labelsize = 10)
    ax.tick_params(axis = "y", labelsize = 10)
    ax.tick_params(axis = "z", labelsize = 10)
    ax.xaxis.get_offset_text().set_fontsize(10)
    ax.yaxis.get_offset_text().set_fontsize(10)
    ax.zaxis.get_offset_text().set_fontsize(10)
    current_figure = Plotter.gcf()

    if save_figure
        Plotter.savefig(title_figure)
    end

    return display(current_figure)

end # main

end # module
