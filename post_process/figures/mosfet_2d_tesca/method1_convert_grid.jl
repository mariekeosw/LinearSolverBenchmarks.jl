module Method1ConvertGrid

using LinearSolverBenchmarks: Mosfet            # to generate the grid
using ExtendableGrids                           # to convert the grid into .dom
using GridVisualize                              # to plot the grid

# Method A: Write grid into TeSCA file format .dom
# How to use this function?
# 1) includet("post_process/figures/mosfet_2d_tesca/method1_convert_grid.jl")
# 2) add Plotter for checking
# 3) Method1ConvertGrid.write_tesca_grid(0, "grid_ref0.dom")
function write_tesca_grid(grid_refinement, filename::String; Plotter = nothing)

    grid = Mosfet.generate_2Dgrid(grid_refinement)

    if Plotter !== nothing
        Plotter.pygui(true) # Plots as Pop-Ups
        GridVisualize.gridplot(grid, Plotter = Plotter, resolution = (600, 450), title = "Grid", xlabel = "Width [m]", ylabel = "Height [m]")
    end

    # convert to "real" \mum
    grid[Coordinates] *= 1.0e6
    trim!(grid)

    # 1 Gate, -2 Drain, -1 Source, -3 Bulk, 999 Neumann
    Base.write(filename, grid, format = "dom", domcodes = (1, -2, -1, -3, 999))

    return grid
end # function tesca

# Method B: Write TeSCA file to simplexgrid format
# How to use this function?
# 1) includet("post_process/figures/mosfet_2d_tesca/method1_convert_grid.jl")
# 2) Method1ConvertGrid.write_ct_grid("grid_ref0.dom")
function write_ct_grid(domfile_tesca::String)

    grid = simplexgrid(domfile_tesca, format = "dom")
    grid[Coordinates] *= 1.0e-6
    trim!(grid)

    return grid
end # function ct

end # module
