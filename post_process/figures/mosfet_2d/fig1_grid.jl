module Fig1Grid

using LinearSolverBenchmarks: Mosfet
using GridVisualize
using LaTeXStrings

# How to use this function?
# 1) includet("post_process/figures/mosfet_2d/fig1_grid.jl")
# 2) using PythonPlot
# 3) Fig1Grid.plot_grid(Plotter = PythonPlot)
function plot_grid(; Plotter = nothing, grid_refinement = 0, save_figure = false, title_figure = nothing)

    if isnothing(Plotter)
        error("Please add PythonPlot.jl in local environment and use Plotter = PythonPlot")
    end

    # generate 2D grid
    grid = Mosfet.generate_2Dgrid(grid_refinement)

    # plot 2D grid
    GridVisualize.gridplot(
        grid, Plotter = Plotter,
        resolution = (800, 500),
        fontsize = 35,
        xlabel = L"\mathrm{width\ [m]}",
        ylabel = L"\mathrm{height\ [m]}"
    )

    # save 2D grid
    if save_figure
        Plotter.savefig(title_figure)
    end

    return nothing
end # function

end # module
