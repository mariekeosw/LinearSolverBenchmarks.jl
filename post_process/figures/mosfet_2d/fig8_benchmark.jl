module Fig8BenchmarkSolver

using JLD2
using Printf
using LaTeXStrings
using ColorSchemes

# How to use this function?
# 1) includet("post_process/figures/mosfet_2d/fig8_benchmark.jl")
# 2) using Plots
# 3) Fig8BenchmarkSolver.plot_benchmark("jld2_files_server/mosfet/mosfet_2d_gmres_bicgstab_19.jld2", Plotter=Plots)
function plot_benchmark(filename; Plotter = nothing, save_figure = false, title_figure = nothing)
    if isnothing(Plotter)
        error("Please add Plots.jl in local environment and use Plotter = Plots")
    end

    # Initialize variables centrally to ensure they are visible in the plotting scope local labels::Vector{String}
    local number_nodes::Dict{String, Any}
    local storage_tlinsolve::Dict{String, Any}

    data = JLD2.load(filename)

    labels = data["labels"]
    number_nodes = data["number_nodes"]
    storage_tlinsolve = data["storage_tlinsolve"]

    # Skip Problemsizes, this is just applied in the very large 2D comparison
    skip_sizes = Set([945576, 757302, 589908, 443394, 317760, 213006, 129132, 66138])

    cols = ColorSchemes.seaborn_deep

    unique_ticks_all = number_nodes[labels[1]]
    keep_tick_idx = findall(N -> !(Int(round(3 * N)) in skip_sizes), unique_ticks_all)
    unique_ticks = unique_ticks_all[keep_tick_idx]
    tick_labels = [string(Int(round(3 * val))) for val in unique_ticks]

    p = Plotter.plot(;
        ylabel = L"\mathrm{runtime\ [s]}",
        xlabel = L"\mathrm{number\ of\ unknowns\ } 3N",
        yscale = :log10,
        xscale = :log10,
        xticks = (unique_ticks, tick_labels),
        #ylims = (0.15e1, 0.75e4), # Figure 6.9
        #ylims = (0.12e1, 1.0e4), # Figure 6.10
        ylims = (0.18e1, 1.8e4), # Figure 6.10
        yticks = 10.0 .^ (0:1:4),
        xrotation = 60,
        title = "",
        legend = :topleft,
        bottom_margin = 10 * Plotter.px,
        guidefontsize = 16,
        framestyle = :box,
        foreground_color_border = :black,
        grid = true,
        gridlinewidth = 0.15,
        gridcolor = :black,
        gridalpha = 0.6,
    )

    for (i, label) in enumerate(labels)

        x = number_nodes[label]
        y = storage_tlinsolve[label]
        keep = findall(N -> !(Int(round(3 * N)) in skip_sizes), x)

        Plotter.plot!(
            p,
            x[keep],
            y[keep];
            linecolor = cols[i],
            linestyle = :solid,
            linewidth = 1.5,
            marker = :circle,
            markersize = 3,
            markerstrokewidth = 0,
            markercolor = cols[i],
            label = label,
        )

    end

    # -------------------------------------------------
    # Complexity / Reference Lines
    # -------------------------------------------------
    N_min = minimum(unique_ticks)
    N_max = maximum(unique_ticks)
    N_ref = collect(range(N_min, N_max, length = 200))

    # One free line in legend
    Plotter.scatter!(p, [], []; label = " ", ms = 0, mc = :white, msc = :white)

    # O(N)
    Plotter.plot!(p, N_ref, 3.0e-2 .* N_ref; color = :black, linestyle = :dash, linewidth = 1, alpha = 0.8, label = L"\mathcal{O}(N)")
    Plotter.plot!(p, N_ref, 1.0e-3 .* N_ref; color = :black, linestyle = :dash, linewidth = 1, alpha = 0.8, primary = false)

    # O(N^2)
    # set new N_ref depending on data
    #N_cut = unique_ticks[round(Int, 0.5 * length(unique_ticks))]
    #N_ref = collect(range(N_min, N_cut, length = 100))
    Plotter.plot!(p, N_ref, 0.8e-5 .* N_ref .^ 2; color = :black, linestyle = :dot, linewidth = 1, alpha = 0.8, label = L"\mathcal{O}(N^2)")

    if save_figure
        Plotter.savefig(p, title_figure)
    end

    return p
end # plot_benchmark

# How to use this function?
# 1) includet("post_process/figures/mosfet_2d/fig8_benchmark.jl")
# 2) using Plots
# 3) Fig8BenchmarkSolver.plot_benchmark_two_files("jld2_files_server/mosfet/mosfet_2d_direct_cn_7.jld2", "jld2_files_server/mosfet/mosfet_2d_direct_le_7.jld2", Plotter=Plots, save_figure=true, title_figure="mosfet_2d_direct_cn_le_7_tlinsolve_2nd.png")
function plot_benchmark_two_files(filename1, filename2; Plotter = nothing, save_figure = false, title_figure = nothing)

    if isnothing(Plotter)
        error("Please add Plots.jl in local environment and use Plotter = Plots")
    end

    # Initialize variables centrally to ensure they are visible in the plotting scope    local labels::Vector{String}
    local number_nodes::Dict{String, Any}
    local storage_tlinsolve::Dict{String, Any}

    # Load both data files into dictionaries
    data1 = JLD2.load(filename1)
    data2 = JLD2.load(filename2)

    # Define the desired sequence and the shortened, cleaner labels - this needs to be changed for new plots
    # Format: (Source_File_ID, Old_Label_In_File, New_Label_For_Plot)
    label_blueprint = [
        (2, "UMFPACK factorization", "UMFPACK (le)"),
        (1, "UMFPACK factorization", "UMFPACK (cn)"),
        (2, "KLU factorization", "KLU (le)"),
        (1, "KLU factorization", "KLU (cn)"),
        (2, "MKLPardiso factorization", "MKLPardiso (le)"),
        (1, "MKLPardiso factorization", "MKLPardiso (cn)"),
        (2, "Sparspak factorization", "Sparspak (le)"),
        (1, "Sparspak factorization", "Sparspak (cn)"),
    ]

    labels = String[]
    number_nodes = Dict{String, Any}()
    storage_tlinsolve = Dict{String, Any}()

    for (file_id, old_label, new_label) in label_blueprint
        push!(labels, new_label)

        if file_id == 1
            number_nodes[new_label] = data1["number_nodes"][old_label]
            storage_tlinsolve[new_label] = data1["storage_tlinsolve"][old_label]
        else
            number_nodes[new_label] = data2["number_nodes"][old_label]
            storage_tlinsolve[new_label] = data2["storage_tlinsolve"][old_label]
        end
    end

    cols = ColorSchemes.seaborn_deep
    unique_ticks = number_nodes[labels[1]]
    tick_labels = [string(Int(round(3 * val))) for val in unique_ticks]

    p = Plotter.plot(;
        ylabel = L"\mathrm{runtime\ [s]}",
        xlabel = L"\mathrm{number\ of\ unknowns\ } 3N",
        yscale = :log10,
        xscale = :log10,
        ylims = (0.25e1, 0.6e4),
        xticks = (unique_ticks, tick_labels),
        title = "",
        legend = :topleft,
        bottom_margin = 10 * Plotter.px,
        guidefontsize = 16,
        xrotation = 60,
        framestyle = :box,
        foreground_color_border = :black,
        linewidth = 2,
        grid = true,
        gridlinewidth = 0.15,
        gridcolor = :black,
        gridalpha = 0.6,
    )

    for (i, label) in enumerate(labels)
        Plotter.plot!(
            p,
            number_nodes[label],
            storage_tlinsolve[label];
            linecolor = cols[i],
            linestyle = :solid,
            marker = :circle,
            markersize = 3,
            markerstrokewidth = 0,
            markercolor = cols[i],
            label = label,
        )

        if i % 2 == 0 && i < length(labels)
            Plotter.scatter!(p, [], []; label = " ", ms = 0, mc = :white, msc = :white)
        end
    end

    # -------------------------------------------------
    # Complexity / Reference Lines
    # -------------------------------------------------
    N_min = minimum(unique_ticks)
    N_max = maximum(unique_ticks)
    N_ref = collect(range(N_min, N_max, length = 200))

    # One free line in legend
    Plotter.scatter!(p, [], []; label = " ", ms = 0, mc = :white, msc = :white)

    # O(N)
    Plotter.plot!(p, N_ref, 3.0e-2 .* N_ref; color = :black, linestyle = :dash, linewidth = 1, alpha = 0.8, label = L"\mathcal{O}(N)")
    Plotter.plot!(p, N_ref, 1.0e-3 .* N_ref; color = :black, linestyle = :dash, linewidth = 1, alpha = 0.8, primary = false)

    # O(N^2)
    Plotter.plot!(p, N_ref, 0.5e-5 .* N_ref .^ 2; color = :black, linestyle = :dot, linewidth = 1, alpha = 0.8, label = L"\mathcal{O}(N^2)")

    if save_figure
        Plotter.savefig(p, title_figure)
    end

    return p
end # plot_benchmark_two_files

end # module
