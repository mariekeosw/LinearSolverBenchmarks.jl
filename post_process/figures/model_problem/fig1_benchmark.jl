module Fig1Benchmark

using JLD2
using Printf
using LaTeXStrings
using ColorSchemes

# How to use these function?
# 0) edit labels + filter algorithms + check factor for reference scaling
# 1) includet("post_process/figures/model_problem/fig1_benchmark.jl")
# 2) using Plots
# 3) Fig1Benchmark.plot_benchmark("jld2_files_server/sym_model_problem/xxx.jld2", Plotter=Plots, save_figure=true, title_figure="xxx")
function plot_benchmark(filename; Plotter = nothing, save_figure = false, title_figure = nothing)

    if isnothing(Plotter)
        error("Please add Plots.jl in local environment and use Plotter = Plots")
    end

    @load filename labels dimensions refinement problem_sizes tests residuals walltimes tlinsolves tlinsolve_setups tlinsolve_solves

    # --- FILTER ALGORITHMS ---
    # Nonsymmetric Krylov Solver: CGS, BiCGSTAB, GMRES
    #aktive_algos = ["MKLPardiso", "CGS-PyAMG-SA", "BiCGSTAB-PyAMG-SA", "GMRES-PyAMG-SA"]
    #aktive_algos = ["MKLPardiso", "BiCGSTAB-PyAMG-RS", "GMRES-PyAMG-RS"]
    #aktive_algos = ["MKLPardiso", "CGS-AMG-SA", "BiCGSTAB-AMG-SA", "GMRES-AMG-SA"]
    #aktive_algos = ["MKLPardiso", "CGS-AMG-RS", "BiCGSTAB-AMG-RS", "GMRES-AMG-RS"]
    #aktive_algos = ["MKLPardiso", "CGS-AMGCL-RS-SPAI0", "BiCGSTAB-AMGCL-RS-SPAI0", "GMRES-AMGCL-RS-SPAI0"]
    #aktive_algos = ["MKLPardiso", "CGS-AMGCL-SA-SPAI0", "BiCGSTAB-AMGCL-SA-SPAI0", "GMRES-AMGCL-SA-SPAI0"]
    ##aktive_algos = ["MKLPardiso", "CGS-AMGCL-RS-ILU0", "BiCGSTAB-AMGCL-RS-ILU0", "GMRES-AMGCL-RS-ILU0"]
    #aktive_algos = ["MKLPardiso", "CGS-AMGCL-SA-ILU0", "BiCGSTAB-AMGCL-SA-ILU0", "GMRES-AMGCL-SA-ILU0"]
    #aktive_algos = ["MKLPardiso", "CGS-ILU0", "BiCGSTAB-ILU0", "GMRES-ILU0"]
    #aktive_algos = ["MKLPardiso", "CGS-ILU-drop-0.1", "BiCGSTAB-ILU-drop-0.1", "GMRES-ILU-drop-0.1"]
    #aktive_algos = ["MKLPardiso", "CGS-ILU-drop-0.01", "BiCGSTAB-ILU-drop-0.01", "GMRES-ILU-drop-0.01"]
    #aktive_algos = ["MKLPardiso", "CGS-LU", "BiCGSTAB-LU", "GMRES-LU"]
    #aktive_algos = ["MKLPardiso", "CGS", "BiCGSTAB", "GMRES"]

    # Symmetric Krylov Solver: CG
    #aktive_algos = ["MKLPardiso", "CG-AMG-RS", "CG-AMG-SA"]
    #aktive_algos = ["MKLPardiso", "CG-AMGCL-SA-ILU0", "CG-AMGCL-RS-ILU0"]
    #aktive_algos = ["MKLPardiso", "CG-AMGCL-SA-SPAI0", "CG-AMGCL-RS-SPAI0"]
    #aktive_algos = ["MKLPardiso", "CG-AMGCL-SA-SPAI0", "CG-AMGCL-RS-SPAI0", "CG-AMGCL-SA-ILU0", "CG-AMGCL-RS-ILU0"]
    #aktive_algos = ["MKLPardiso", "CG-ILU0", "CG-ILU-drop-0.01", "CG-ILU-drop-0.1"]
    #aktive_algos = ["MKLPardiso", "CG-LU"]
    #aktive_algos = ["MKLPardiso", "CG"]

    #labels = filter(l -> l in aktive_algos, labels)

    cols = ColorSchemes.seaborn_deep
    lst = [:dot, :dash, :solid]

    p = Plotter.plot(;
        ylabel = L"\mathrm{runtime\ [s]}",
        xlabel = L"\mathrm{problem\ size}",
        yscale = :log10,
        xscale = :log10,
        ylims = (1.0e-5, 1.0e3), # (1.0e-5, 0.2e1) for 2D
        #xlims = ,
        yticks = 10.0 .^ (-5:1:3),
        xticks = 10.0 .^ (1:8),
        bottom_margin = 10 * Plotter.px,
        legend = :topleft,
        guidefontsize = 16,
        framestyle = :box,
        foreground_color_border = :black,
        linewidth = 2,
        grid = true,
        gridlinewidth = 0.15,
        gridcolor = :black,
        gridalpha = 0.2,
    )

    dimensions = [3]

    for dimension in dimensions
        for (i, label) in enumerate(labels)
            key = (dimension, label)

            valid = (problem_sizes[key] .> 0) .& isfinite.(tlinsolves[key])    # needed for GMRES PyAMG (large problem sizes did not run through and consequently have problemsize=0)

            Plotter.plot!(
                p,
                problem_sizes[key][valid],
                tlinsolves[key][valid];
                linecolor = cols[i],
                linestyle = lst[dimension],
                marker = :circle,
                markersize = 3,
                markerstrokewidth = 0,
                markercolor = cols[i],
                primary = false,
            )
        end
    end

    # -------------------------------------------------
    # Custom Legend
    # -------------------------------------------------
    # Custom legend: colors represent algorithms
    for (i, label) in enumerate(labels)

        Plotter.plot!(
            p,
            Float64[],
            Float64[];
            color = cols[i],
            linestyle = :solid,
            linewidth = 1.5,
            marker = :circle,
            markersize = 3,
            markerstrokewidth = 0,
            label = label,
        )
    end

    # One free line in legend
    Plotter.scatter!(p, [], []; label = " ", ms = 0, mc = :white, msc = :white)

    # Custom legend: line styles represent dimensions
    for dimension in dimensions
        Plotter.plot!(
            p,
            Float64[],
            Float64[];
            color = :black,
            linestyle = lst[dimension],
            linewidth = 1.5,
            label = "$(dimension)D",
        )
    end

    # -------------------------------------------------
    # Complexity / Reference Lines
    # -------------------------------------------------
    # Collect all problem sizes
    all_N = vcat(values(problem_sizes)...)
    all_N = all_N[all_N .> 0] # we need this in Consequently,case GMRES PyAMG

    # Dynamic range
    N_min = minimum(all_N)
    N_max = maximum(all_N)

    # Smooth reference grid
    N_ref = collect(range(N_min, N_max, length = 200))

    # One free line in legend
    Plotter.scatter!(p, [], []; label = " ", ms = 0, mc = :white, msc = :white)

    # O(N)
    Plotter.plot!(p, N_ref, 1.0e-4 .* N_ref; color = :black, linestyle = :dash, linewidth = 1, alpha = 0.8, label = L"\mathcal{O}(N)")
    Plotter.plot!(p, N_ref, 1.0e-7 .* N_ref; color = :black, linestyle = :dash, linewidth = 1, alpha = 0.8, primary = false)

    # O(N^2)
    Plotter.plot!(p, N_ref, 1.0e-8 .* N_ref .^ 2; color = :black, linestyle = :dot, linewidth = 1, alpha = 0.8, label = L"\mathcal{O}(N^2)")

    if save_figure
        Plotter.savefig(p, title_figure)
    end

    return p
end # plot one file

function plot_benchmark_two_files(
        filename1, filename2;
        Plotter = nothing,
        save_figure = false,
        title_figure = nothing,
        active_dimensions = nothing
    )

    if isnothing(Plotter)
        error("Please add Plots.jl in local environment and use Plotter = Plots")
    end

    data1 = JLD2.load(filename1)
    data2 = JLD2.load(filename2)

    # Combine labels and dimensions from both files
    labels = unique(vcat(data1["labels"], data2["labels"]))

    dimensions = sort(unique(vcat(data1["dimensions"], data2["dimensions"])))

    # --- FILTER DIMENSIONS ---
    if !isnothing(active_dimensions)
        selected_dimensions = active_dimensions isa Integer ?
            [active_dimensions] :
            active_dimensions

        dimensions = filter(
            dimension -> dimension in selected_dimensions,
            dimensions,
        )

        isempty(dimensions) &&
            error("None of the selected dimensions occur in the benchmark file.")
    end

    # Combine dictionaries containing the actual benchmark results
    problem_sizes = merge(data1["problem_sizes"], data2["problem_sizes"])

    tlinsolves = merge(data1["tlinsolves"], data2["tlinsolves"])

    # --- FILTER ALGORITHMS ---
    # Direkt Solver
    #aktive_algos = ["MKLPardiso", "KLU", "UMFPACK", "Sparspak"]

    # Nonsymmetric Krylov Solver: CGS, BiCGSTAB, GMRES
    aktive_algos = ["MKLPardiso", "BiCGSTAB-PyAMG-SA", "GMRES-PyAMG-SA"]
    #aktive_algos = ["MKLPardiso", "BiCGSTAB-PyAMG-RS", "GMRES-PyAMG-RS"]
    #aktive_algos = ["MKLPardiso", "CGS-AMG-SA", "BiCGSTAB-AMG-SA", "GMRES-AMG-SA"]
    #aktive_algos = ["MKLPardiso", "CGS-AMG-RS", "BiCGSTAB-AMG-RS", "GMRES-AMG-RS"]
    #aktive_algos = ["MKLPardiso", "CGS-AMGCL-RS-SPAI0", "BiCGSTAB-AMGCL-RS-SPAI0", "GMRES-AMGCL-RS-SPAI0"]
    #aktive_algos = ["MKLPardiso", "CGS-AMGCL-SA-SPAI0", "BiCGSTAB-AMGCL-SA-SPAI0", "GMRES-AMGCL-SA-SPAI0"]
    #aktive_algos = ["MKLPardiso", "CGS-AMGCL-RS-ILU0", "BiCGSTAB-AMGCL-RS-ILU0", "GMRES-AMGCL-RS-ILU0"]
    #aktive_algos = ["MKLPardiso", "CGS-AMGCL-SA-ILU0", "BiCGSTAB-AMGCL-SA-ILU0", "GMRES-AMGCL-SA-ILU0"]
    #aktive_algos = ["MKLPardiso", "CGS-ILU0", "BiCGSTAB-ILU0", "GMRES-ILU0"]
    #aktive_algos = ["MKLPardiso", "CGS-ILU-drop-0.1", "BiCGSTAB-ILU-drop-0.1", "GMRES-ILU-drop-0.1"]
    #aktive_algos = ["MKLPardiso", "CGS-ILU-drop-0.01", "BiCGSTAB-ILU-drop-0.01", "GMRES-ILU-drop-0.01"]
    #aktive_algos = ["MKLPardiso", "CGS", "BiCGSTAB", "GMRES"]

    # Symmetric Krylov Solver: CG
    #aktive_algos = ["CG", "CG-PyAMG-RS", "CG-ILU0"]
    #aktive_algos = ["CG", "CG-AMG-RS", "CG-ILU0"]
    #aktive_algos = ["CG", "CG-AMGCL-SA-SPAI0", "CG-ILU0"]
    #aktive_algos = ["MKLPardiso", "CG-PyAMG-RS", "CG-PyAMG-SA"]
    #aktive_algos = ["MKLPardiso", "CG-AMG-RS", "CG-AMG-SA"]
    #aktive_algos = ["MKLPardiso", "CG-AMGCL-SA-ILU0", "CG-AMGCL-RS-ILU0"]
    #aktive_algos = ["MKLPardiso", "CG-AMGCL-SA-SPAI0", "CG-AMGCL-RS-SPAI0"]
    #aktive_algos = ["MKLPardiso", "CG-AMGCL-SA-SPAI0", "CG-AMGCL-RS-SPAI0", "CG-AMGCL-SA-ILU0", "CG-AMGCL-RS-ILU0"]
    #aktive_algos = ["MKLPardiso", "CG-ILU0", "CG-ILU-drop-0.01", "CG-ILU-drop-0.1"]
    #aktive_algos = ["MKLPardiso", "CG"]

    labels = filter(l -> l in aktive_algos, labels)

    cols = ColorSchemes.seaborn_deep
    lst = [:dot, :dash, :solid]
    #lst = [:solid, :solid, :solid] # if just one dimension considered

    p = Plotter.plot(;
        ylabel = L"\mathrm{runtime\ [s]}",
        xlabel = L"\mathrm{number\ of\ unknowns\ } N",
        yscale = :log10,
        xscale = :log10,
        ylims = (1.0e-4, 1.0e4),
        xlims = (1.0e2 / 1.3, 42144192 * 1.3), # ref 22 max N ist : 42144192 # ref 16 max N ist: 328509
        yticks = 10.0 .^ (-4:2:4),
        xticks = 10.0 .^ (2:8),
        left_margin = 6 * Plotter.px,
        bottom_margin = 3 * Plotter.px,
        legend = :topleft,
        guidefontsize = 16,
        framestyle = :box,
        foreground_color_border = :black,
        linewidth = 2,
        grid = true,
        gridlinewidth = 0.15,
        gridcolor = :black,
        gridalpha = 0.6,
    )

    for dimension in dimensions
        for (i, label) in enumerate(labels)
            key = (dimension, label)

            valid = (problem_sizes[key] .> 0) .& isfinite.(tlinsolves[key])    # needed for GMRES PyAMG (large problem sizes did not run through and consequently have problemsize=0)

            Plotter.plot!(
                p,
                problem_sizes[key][valid],
                tlinsolves[key][valid];
                linecolor = cols[i],
                linestyle = lst[dimension],
                marker = :circle,
                markersize = 3,
                markerstrokewidth = 0,
                markercolor = cols[i],
                primary = false,
            )
        end
    end

    # -------------------------------------------------
    # Custom Legend
    # -------------------------------------------------
    # Custom legend: colors represent algorithms
    for (i, label) in enumerate(labels)

        Plotter.plot!(
            p,
            Float64[],
            Float64[];
            color = cols[i],
            linestyle = :solid,
            linewidth = 1.5,
            marker = :circle,
            markersize = 3,
            markerstrokewidth = 0,
            label = label,
        )
    end

    # One free line in legend
    Plotter.scatter!(p, [], []; label = " ", ms = 0, mc = :white, msc = :white)

    # Custom legend: line styles represent dimensions
    for dimension in dimensions
        Plotter.plot!(
            p,
            Float64[],
            Float64[];
            color = :black,
            linestyle = lst[dimension],
            linewidth = 1.5,
            label = "$(dimension)D",
        )
    end

    # -------------------------------------------------
    # Complexity / Reference Lines
    # -------------------------------------------------
    # Collect all problem sizes
    all_N = vcat(values(problem_sizes)...)
    all_N = all_N[all_N .> 0] # we need this in case GMRES PyAMG

    # Dynamic range
    N_min = minimum(all_N)
    N_max = 42144192 # maximum(all_N) # ref 22 max N ist : 42144192 # ref 16 max N ist: 328509

    # Smooth reference grid
    N_ref = collect(range(N_min, N_max, length = 200))

    # One free line in legend
    Plotter.scatter!(p, [], []; label = " ", ms = 0, mc = :white, msc = :white)

    # All Dimensions + Direct Solver
    # O(N)
    Plotter.plot!(p, N_ref, 1.0e-4 .* N_ref; color = :black, linestyle = :dash, linewidth = 1, alpha = 0.8, label = L"\mathcal{O}(N)")
    Plotter.plot!(p, N_ref, 1.0e-7 .* N_ref; color = :black, linestyle = :dash, linewidth = 1, alpha = 0.8, primary = false)
    # O(N^2)
    # Change for the direct plot
    #Plotter.plot!(p, N_ref, 0.5e-7 .* N_ref .^ 2; color = :black, linestyle = :dot, linewidth = 1, alpha = 0.8, label = L"\mathcal{O}(N^2)")
    Plotter.plot!(p, N_ref, 1.0e-9 .* N_ref .^ 2; color = :black, linestyle = :dot, linewidth = 1, alpha = 0.8, label = L"\mathcal{O}(N^2)")

    if save_figure
        Plotter.savefig(p, title_figure)
    end

    return p
end # plot two files

function plot_benchmark_three_files(
        filename1, filename2, filename3;
        Plotter = nothing,
        save_figure = false,
        title_figure = nothing,
        active_dimensions = nothing
    )

    if isnothing(Plotter)
        error("Please add Plots.jl in local environment and use Plotter = Plots")
    end

    data1 = JLD2.load(filename1)
    data2 = JLD2.load(filename2)
    data3 = JLD2.load(filename3)

    # Combine labels and dimensions from both files
    labels = unique(vcat(data1["labels"], data2["labels"], data3["labels"]))

    dimensions = sort(unique(vcat(data1["dimensions"], data2["dimensions"], data3["dimensions"])))

    # --- FILTER DIMENSIONS ---
    if !isnothing(active_dimensions)
        selected_dimensions = active_dimensions isa Integer ?
            [active_dimensions] :
            active_dimensions

        dimensions = filter(
            dimension -> dimension in selected_dimensions,
            dimensions,
        )

        isempty(dimensions) &&
            error("None of the selected dimensions occur in the benchmark file.")
    end

    # Combine dictionaries containing the actual benchmark results
    problem_sizes = merge(data1["problem_sizes"], data2["problem_sizes"], data3["problem_sizes"])

    tlinsolves = merge(data1["tlinsolves"], data2["tlinsolves"], data3["tlinsolves"])

    # --- FILTER ALGORITHMS ---
    # Nonsymmetric Krylov Solver: CGS, BiCGSTAB, GMRES
    #aktive_algos = ["MKLPardiso", "CGS", "BiCGSTAB", "GMRES"]
    aktive_algos = ["MKLPardiso", "CGS-AMG-SA", "BiCGSTAB-AMG-SA", "GMRES-AMG-SA"]
    #aktive_algos = ["MKLPardiso", "CGS-AMG-RS", "BiCGSTAB-AMG-RS", "GMRES-AMG-RS"]
    #aktive_algos = ["MKLPardiso", "CGS-AMGCL-SA-ILU0", "BiCGSTAB-AMGCL-SA-ILU0", "GMRES-AMGCL-SA-ILU0"]
    #aktive_algos = ["MKLPardiso", "CGS-AMGCL-RS-ILU0", "BiCGSTAB-AMGCL-RS-ILU0", "GMRES-AMGCL-RS-ILU0"]

    # Symmetric Krylov Solver: CG
    #aktive_algos = ["MKLPardiso", "CG-AMG-RS", "CG-AMG-SA"]
    #aktive_algos = ["CG", "CG-PyAMG-RS", "CG-ILU0"]
    #aktive_algos = ["CG", "CG-AMG-RS", "CG-ILU0"]
    #aktive_algos = ["CG", "CG-AMGCL-SA-SPAI0", "CG-ILU0"]

    labels = filter(l -> l in aktive_algos, labels)

    cols = ColorSchemes.seaborn_deep
    #lst = [:solid, :solid, :solid]
    lst = [:dot, :dash, :solid]

    p = Plotter.plot(;
        ylabel = L"\mathrm{runtime\ [s]}",
        xlabel = L"\mathrm{number\ of\ unknowns\ } N",
        yscale = :log10,
        xscale = :log10,
        ylims = (1.0e-4, 1.0e4),
        xlims = (1.0e2 / 1.3, 42144192 * 1.3), # ref 22 N_max: 42144192
        yticks = 10.0 .^ (-4:2:4),
        xticks = 10.0 .^ (2:8),
        left_margin = 6 * Plotter.px,
        bottom_margin = 3 * Plotter.px,
        legend = :topleft,
        guidefontsize = 16,
        framestyle = :box,
        foreground_color_border = :black,
        linewidth = 2,
        grid = true,
        gridlinewidth = 0.15,
        gridcolor = :black,
        gridalpha = 0.6,
    )

    for dimension in dimensions
        for (i, label) in enumerate(labels)
            key = (dimension, label)

            valid = (problem_sizes[key] .> 0) .& isfinite.(tlinsolves[key])    # needed for GMRES PyAMG (large problem sizes did not run through and consequently have problemsize=0)

            Plotter.plot!(
                p,
                problem_sizes[key][valid],
                tlinsolves[key][valid];
                linecolor = cols[i],
                linestyle = lst[dimension],
                marker = :circle,
                markersize = 3,
                markerstrokewidth = 0,
                markercolor = cols[i],
                primary = false,
            )
        end
    end

    # -------------------------------------------------
    # Custom Legend
    # -------------------------------------------------
    # Custom legend: colors represent algorithms
    for (i, label) in enumerate(labels)

        Plotter.plot!(
            p,
            Float64[],
            Float64[];
            color = cols[i],
            linestyle = :solid,
            linewidth = 1.5,
            marker = :circle,
            markersize = 3,
            markerstrokewidth = 0,
            label = label,
        )
    end

    # One free line in legend
    Plotter.scatter!(p, [], []; label = " ", ms = 0, mc = :white, msc = :white)

    # Custom legend: line styles represent dimensions
    for dimension in dimensions
        Plotter.plot!(
            p,
            Float64[],
            Float64[];
            color = :black,
            linestyle = lst[dimension],
            linewidth = 1.5,
            label = "$(dimension)D",
        )
    end

    # -------------------------------------------------
    # Complexity / Reference Lines
    # -------------------------------------------------
    # Collect all problem sizes
    all_N = vcat(values(problem_sizes)...)
    all_N = all_N[all_N .> 0] # we need this in case GMRES PyAMG

    # Dynamic range
    N_min = minimum(all_N)
    N_max = 42144192 #maximum(all_N)

    # Smooth reference grid
    N_ref = collect(range(N_min, N_max, length = 200))

    # One free line in legend
    Plotter.scatter!(p, [], []; label = " ", ms = 0, mc = :white, msc = :white)

    # 2D
    # O(N)
    #Plotter.plot!(p, N_ref, 1.0e-6 .* N_ref; color = :black, linestyle = :dash, linewidth = 1, alpha = 0.8, label = L"\mathcal{O}(N)")
    # O(N^(3/2))
    #Plotter.plot!(p, N_ref, 1.0e-8 .* N_ref .^ (3 / 2); color = :black, linestyle = :dashdotdot, linewidth = 1, alpha = 0.8, label = L"\mathcal{O}(N^{3/2})")

    # 3D
    # O(N)
    #Plotter.plot!(p, N_ref, 1.0e-6 .* N_ref; color = :black, linestyle = :dash, linewidth = 1, alpha = 0.8, label = L"\mathcal{O}(N)")
    # O(N^(4/3))
    #Plotter.plot!(p, N_ref, 1.0e-7 .* N_ref .^ (4 / 3); color = :black, linestyle = :dashdot, linewidth = 1, alpha = 0.8, label = L"\mathcal{O}(N^{4/3})")

    # All Dimensions + Direct Solver
    # O(N)
    Plotter.plot!(p, N_ref, 1.0e-4 .* N_ref; color = :black, linestyle = :dash, linewidth = 1, alpha = 0.8, label = L"\mathcal{O}(N)")
    Plotter.plot!(p, N_ref, 1.0e-7 .* N_ref; color = :black, linestyle = :dash, linewidth = 1, alpha = 0.8, primary = false)
    # O(N^2)
    # Change for the direct plot
    #Plotter.plot!(p, N_ref, 0.5e-7 .* N_ref .^ 2; color = :black, linestyle = :dot, linewidth = 1, alpha = 0.8, label = L"\mathcal{O}(N^2)")
    Plotter.plot!(p, N_ref, 1.0e-9 .* N_ref .^ 2; color = :black, linestyle = :dot, linewidth = 1, alpha = 0.8, label = L"\mathcal{O}(N^2)")

    if save_figure
        Plotter.savefig(p, title_figure)
    end

    return p
end # plot three files

end # module
