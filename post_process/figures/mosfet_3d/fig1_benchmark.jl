module Fig1BenchmarkSolver

using JLD2
using Printf
using LaTeXStrings
using ColorSchemes

# How to use this function?
# 1) includet("post_process/figures/mosfet_3d/fig1_benchmark.jl")
# 2) using Plots
# 3) Fig1BenchmarkSolver.plot_benchmark_two_files("jld2_files_server/mosfet/mosfet_3d_pardiso_1.jld2", "jld2_files_server/mosfet/mosfet_3d_gmres_bicgstab_2.jld2", Plotter = Plots)
function plot_benchmark_two_files(filename1, filename2; Plotter = nothing, save_figure = false, title_figure = nothing)
    if isnothing(Plotter)
        error("Please add Plots.jl in local environment and use Plotter = Plots")
    end

    # Initialize variables centrally to ensure they are visible in the plotting scope local labels::Vector{String}
    local number_nodes::Dict{String, Any}
    local storage_tlinsolve::Dict{String, Any}

    data1 = JLD2.load(filename1)
    data2 = JLD2.load(filename2)

    labels = vcat(data1["labels"], data2["labels"])
    number_nodes = merge(data1["number_nodes"], data2["number_nodes"])
    storage_tlinsolve = merge(data1["storage_tlinsolve"], data2["storage_tlinsolve"])

    cols = ColorSchemes.seaborn_deep

    unique_ticks = number_nodes[labels[2]]
    tick_labels = [string(Int(round(3 * val))) for val in unique_ticks]

    p = Plotter.plot(;
        ylabel = L"\mathrm{runtime\ [s]}",
        xlabel = L"\mathrm{number\ of\ unknowns\ } 3N",
        yscale = :log10,
        xscale = :log10,
        xticks = (unique_ticks, tick_labels),
        yticks = 10.0 .^ (2:1:5),
        ylims = (1.3e2, 1.0e5),
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

        Plotter.plot!(
            p,
            number_nodes[label],
            storage_tlinsolve[label];
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
    Plotter.plot!(p, N_ref, 1.5e-2 .* N_ref; color = :black, linestyle = :dash, linewidth = 1, alpha = 0.8, label = L"\mathcal{O}(N)")
    Plotter.plot!(p, N_ref, 3.0e-3 .* N_ref; color = :black, linestyle = :dash, linewidth = 1, alpha = 0.8, primary = false)
    # O(N^2)
    Plotter.plot!(p, N_ref, 0.4e-6 .* N_ref .^ 2; color = :black, linestyle = :dot, linewidth = 1, alpha = 0.8, label = L"\mathcal{O}(N^2)")

    if save_figure
        Plotter.savefig(p, title_figure)
    end

    return p
end # plot_benchmark

end # module
