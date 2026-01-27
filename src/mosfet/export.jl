################################################################################
# Helperfunction for Plotting
################################################################################
# Took from Dilara: TexturedPerovskiteSolarCells.jl
# THX at JF!!!
# https://github.com/j-fu/GridVisualize.jl/blob/1f2b299a436b7750702ccca282fa14152d80ebf9/src/pyplot.jl#L86
function tridata(grid::ExtendableGrid)
    coord = grid[Coordinates]
    cellnodes = Matrix(grid[CellNodes])
    return coord[1, :], coord[2, :], transpose(cellnodes .- 1)
end

################################################################################
# Heatmap Matrix (with Plots.jl)
################################################################################
# Helperfunction for heatmap
function specieswise_permutation(N)
    return vcat(
        1:3:N,      # φ_n;1, ..., φ_n;N
        2:3:N,      # φ_p;1, ..., φ_p;N
        3:3:N       # ψ_1, ..., ψ_N
    )
end

function plot_heatmap_matrix(system, solution; blocksize = 3, permutation = false, Plotter = nothing)
    # Jacobian
    residual, jacobian_sparse = VoronoiFVM.evaluate_residual_and_jacobian(system, solution)
    jacobian = Matrix(jacobian_sparse)
    N = size(jacobian, 1) # Dimension jacobian

    # Permutation to specieswise ordering
    if permutation
        perm = specieswise_permutation(N)
        jacobian = jacobian[perm, perm]
    end

    # Reduced jacobian (one block -> one entry = max. entry of the block)
    n = cld(N, blocksize) # Dimension reduced jacobian
    jacobian_blockmax = fill(NaN, n, n)

    for i in 1:n
        for j in 1:n
            rows = ((i - 1) * blocksize + 1):min(i * blocksize, N)
            cols = ((j - 1) * blocksize + 1):min(j * blocksize, N)

            block = @view jacobian[rows, cols]
            jacobian_blockmax[i, j] = maximum(abs.(block))
        end
    end

    # Log scale for visualization
    jacobian_blockmax_log = log10.(abs.(jacobian_blockmax) .+ 1.0e-16)
    jacobian_blockmax_log[jacobian_blockmax_log .< -15] .= NaN # Hide very small values

    # Axis scaling to show the actual dimension of jacobian
    tickpos = 20:20:n
    ticklab = string.((tickpos .- 1) .* blocksize .+ 1)

    # Create plot
    plt = Plotter.heatmap(
        jacobian_blockmax_log,
        yflip = true,
        color = Plotter.cgrad(:inferno, rev = true),
        xticks = (tickpos, ticklab),
        yticks = (tickpos, ticklab),
        xlabel = L"\mathrm{column\ index\ } j",
        ylabel = L"\mathrm{row\ index\ } i",
        colorbar_title = L"\log_{10}(|J_{ij}|)",
        colorbar_titlefontsize = 16,
        right_margin = 15 * Plotter.px,
        # Bigger axis labels
        guidefontsize = 16,
        # Outer frame thicker
        framestyle = :box,
        foreground_color_border = :black,
        linewidth = 2,
        # Grid lines
        grid = true,
        gridlinewidth = 0.15,
        gridcolor = :black,
        gridalpha = 0.6,
    )

    display(plt)

    return plt
end

################################################################################
# Latex Export
################################################################################
# Export function to create detailed latex table for all parameter values and runtimes
function export_detailed_benchmark_to_latex(benchmark_data::BenchmarkSeries; filename_out::AbstractString = "benchmark_detailed.tex")
    # Helper: format parameter values
    function format_value(x::Float64; mode::Symbol)
        if mode === :equilibrium
            # compact scientific / plain formatting for λ1
            if x == 0.0
                return "0"
            elseif abs(x) < 1.0e-2 || abs(x) >= 1.0e2
                s = @sprintf("%.0e", x)   # e.g. 1e-20
                base, exp = split(s, "e")
                exp_int = parse(Int, exp)
                if base == "1"
                    return "\$10^{$(exp_int)}\$"
                else
                    return "\$$(base) \\cdot 10^{$(exp_int)}\$"
                end
            else
                return string(x)
            end
        elseif mode === :bias
            # cleaner formatting for gate/drain voltages
            return @sprintf("%.2f", x)
        else
            return string(x)
        end
    end

    # Helper: format runtime
    function format_time(x::Float64)
        if x == 0.0
            return ""
        elseif x < 1.0e-3
            return @sprintf("%.2e", x)
        elseif x < 1.0
            return @sprintf("%.4f", x)
        else
            return @sprintf("%.3f", x)
        end
    end

    # Helper: one table cell triplet
    function format_row_entry(entry::BenchmarkStep, mode::Symbol)
        return (
            format_value(entry.value; mode = mode),
            string(entry.newton_iters),
            format_time(entry.tlinsolve),
        )
    end

    # Helper: empty triplet
    empty_triplet() = ("", "", "")

    nrows = maximum(
        (
            length(benchmark_data.equilibrium),
            length(benchmark_data.gate),
            length(benchmark_data.drain),
        )
    )

    open(filename_out, "w") do io
        println(io, "\\begin{table}[h]")
        println(io, "\\centering")
        println(io, "\\caption{<Enter caption here>}")
        println(io, "\\label{<enter_label_here>}")
        println(io, "\\renewcommand{\\arraystretch}{1.2}")
        println(io)
        println(io, "\\begin{tabular}{ccc@{\\hspace{1.2cm}}ccc@{\\hspace{1.2cm}}ccc}")
        println(io, "\\toprule")
        println(io, "\\multicolumn{3}{c@{\\hspace{1.2cm}}}{\\textbf{In Equilibrium}} & \\multicolumn{3}{c@{\\hspace{1.2cm}}}{\\textbf{Gate-Sweep}} & \\multicolumn{3}{c}{\\textbf{Drain-Sweep}} \\\\")
        println(io, "\\midrule")
        println(io, "\$\\lambda_1\$ & \$n\$ & \$t_{\\text{lin}} [s]\$ & \$V_G\$ & \$n\$ & \$t_{\\text{lin}} [s]\$ & \$V_D\$ & \$n\$ & \$t_{\\text{lin}} [s]\$\\\\")
        println(io, "\\midrule")

        for i in 1:nrows
            eq = i <= length(benchmark_data.equilibrium) ? format_row_entry(benchmark_data.equilibrium[i], :equilibrium) : empty_triplet()
            ga = i <= length(benchmark_data.gate) ? format_row_entry(benchmark_data.gate[i], :bias) : empty_triplet()
            dr = i <= length(benchmark_data.drain) ? format_row_entry(benchmark_data.drain[i], :bias) : empty_triplet()

            println(
                io,
                eq[1], " & ", eq[2], " & ", eq[3], " & ",
                ga[1], " & ", ga[2], " & ", ga[3], " & ",
                dr[1], " & ", dr[2], " & ", dr[3], "\\\\"
            )
        end

        println(io, "\\bottomrule")
        println(io, "\\end{tabular}")
        println(io, "\\end{table}")
    end

    return nothing
end

################################################################################
# Doping Profile
################################################################################
# Region Doping - Constant
function plot_region_doping_2D(ctsys, Plotter)

    grid = ctsys.fvmsys.grid
    data = ctsys.fvmsys.physics.data

    params = data.params

    cellregions = grid[CellRegions]
    cellValue = zeros(length(cellregions))

    for i in eachindex(cellregions)
        # determine doping value in cell and number of cell nodes
        cellValue[i] = (params.doping[1, cellregions[i]] - params.doping[2, cellregions[i]]) * 1.0e-6
    end

    Plotter.figure()
    m = Plotter.tripcolor(
        tridata(grid)...,
        cellValue;
        shading = "flat",
        alpha = 0.75,
        rasterized = true
    )
    # instead of color bar
    vals = sort(unique(cellValue))
    cmap = m.get_cmap()
    norm = m.norm
    Patch = Plotter.matplotlib.patches.Patch

    handles = [Patch(facecolor = cmap(norm(v)), edgecolor = "k") for v in vals]
    labels = ["$(v)" for v in vals]

    Plotter.xlabel("Width [m]", fontsize = 12)
    Plotter.ylabel("Height [m]", fontsize = 12)
    Plotter.title("Doping Profile", fontsize = 16)

    ax = Plotter.gca()
    ax.xaxis.set_major_locator(Plotter.matplotlib.ticker.MultipleLocator(0.5e-6))
    ax.yaxis.set_major_locator(Plotter.matplotlib.ticker.MultipleLocator(0.5e-6))
    ax.tick_params(axis = "both", which = "major", labelsize = 10)
    ax.xaxis.get_offset_text().set_fontsize(10)
    ax.yaxis.get_offset_text().set_fontsize(10)
    ax.legend(handles, labels, title = "Doping values [cm⁻³]", loc = "lower right", fontsize = 10, title_fontsize = 8)

    Plotter.tight_layout()
    current_figure = Plotter.gcf()
    return display(current_figure)
end

# Nodal Doping
function plot_nodal_doping_2D(ctsys, Plotter)

    grid = ctsys.fvmsys.grid

    data = ctsys.fvmsys.physics.data
    paramsnodal = data.paramsnodal
    doping_scaled = (-1) .* paramsnodal.doping .* 1.0e-6

    Plotter.figure()
    Plotter.tripcolor(
        tridata(grid)..., doping_scaled,
        norm = Plotter.matplotlib.colors.SymLogNorm(linthresh = 1.0e12, vmin = -1.0e16, vmax = 1.0e20),
        shading = "gouraud",
        rasterized = true
    )
    Plotter.colorbar(orientation = "vertical", extend = "both")
    Plotter.title("Nodal Doping Profile", fontsize = 16)
    Plotter.xlabel("Width [m]", fontsize = 12)
    Plotter.ylabel("Height [m]", fontsize = 12)
    current_figure = Plotter.gcf()
    return display(current_figure)
end

function plotsurface_nodal_doping_2D(ctsys, Plotter)

    grid = ctsys.fvmsys.grid
    XX = grid[Coordinates][1, :]; YY = grid[Coordinates][2, :]

    data = ctsys.fvmsys.physics.data
    paramsnodal = data.paramsnodal
    doping_scaled = (-1) .* paramsnodal.doping .* 1.0e-6

    fig = Plotter.figure()
    ax = fig.add_subplot(111, projection = "3d")
    ax.plot_trisurf(XX[:], YY[:], doping_scaled)

    Plotter.title("Nodal Doping Profile", fontsize = 16)
    Plotter.xlabel("Width [m]", fontsize = 12)
    Plotter.ylabel("Height [m]", fontsize = 12)
    current_figure = Plotter.gcf()
    return display(current_figure)
end
