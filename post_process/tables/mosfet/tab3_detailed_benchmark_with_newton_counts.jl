module Tab3DetailedBenchmark

using JLD2
using Printf

# How to use this function?
# 1) includet("post_process/tables/mosfet/tab3_detailed_benchmark_with_newton_counts.jl")
# 2) Tab3DetailedBenchmark.export_benchmark_to_latex("jld2_files_server/mosfet/mosfet_3d_gmres_bicgstab_2.jld2", "xxx.tex")
function export_benchmark_to_latex(filename_in, filename_out)

    @load filename_in labels ref number_nodes walltimes storage_tlinsolve storage_tlinsolve_setup storage_tlinsolve_solve storage_newton_iters tests tests_phy

    open(filename_out, "w") do io
        println(io, "% Requires \\usepackage{booktabs}")
        println(io, "% Requires \\usepackage{longtable}")
        println(io)
        println(io, "\\setlength{\\tabcolsep}{3pt}")
        println(io, "\\footnotesize")
        println(io)
        println(io, "\\begin{longtable}{llrrrrrrr}")
        println(io, "\\caption{Detailed benchmark results for the 2D or 3D MOSFET simulations. The Jacobian size is three times the number of grid nodes due to the three coupled unknowns.}")
        println(io, "\\label{tab:benchmark-full}\\\\")
        println(io)

        println(io, "\\toprule")
        println(io, "Method & \$3N\$ & Wall [s] & Linear [s] & Setup [s] & Solve [s] & \$n\$ & Test 1 & Test 2 \\\\")
        println(io, "\\midrule")
        println(io, "\\endfirsthead")
        println(io)

        println(io, "\\toprule")
        println(io, "Method & \$3N\$ & Wall [s] & Linear [s] & Setup [s] & Solve [s] & \$n\$ & Test 1 & Test 2 \\\\")
        println(io, "\\midrule")
        println(io, "\\endhead")
        println(io)

        println(io, "\\midrule")
        println(io, "\\multicolumn{9}{r}{Continued on next page}")
        println(io, "\\endfoot")
        println(io)

        println(io, "\\bottomrule")
        println(io, "\\endlastfoot")
        println(io)

        for i in eachindex(ref)
            n_nodes = number_nodes[labels[1]][i]
            size_jac = 3 * n_nodes

            for label in labels
                wall = walltimes[label][i]
                total = storage_tlinsolve[label][i]
                setup = storage_tlinsolve_setup[label][i]
                solve = storage_tlinsolve_solve[label][i]
                niter = storage_newton_iters[label][i]
                test = tests[label][i]
                test_phy = tests_phy[label][i]

                method = replace(label, "Block Jacobi " => "")
                method = latex_escape(method)

                println(
                    io,
                    method, " & ",
                    size_jac, " & ",
                    @sprintf("%.2f", wall), " & ",
                    @sprintf("%.2f", total), " & ",
                    @sprintf("%.2f", setup), " & ",
                    @sprintf("%.2f", solve), " & ",
                    niter, " & ",
                    @sprintf("%.3f", test), " & ",
                    @sprintf("%.3f", test_phy), " \\\\"
                    #format_sci_latex(test_phy), " \\\\"
                )
            end

            if i != lastindex(ref)
                println(io, "\\\\[0.3em]")
                println(io)
            end
        end

        println(io, "\\end{longtable}")
        println(io, "\\normalsize")
    end

    return nothing
end # function

# Helper function to escape special LaTeX characters in method labels
function latex_escape(s::AbstractString)
    s = replace(s, "\\" => "\\textbackslash{}")
    s = replace(s, "_" => "\\_")
    s = replace(s, "&" => "\\&")
    s = replace(s, "%" => "\\%")
    s = replace(s, "#" => "\\#")
    s = replace(s, "{" => "\\{")
    s = replace(s, "}" => "\\}")
    return s
end

function format_sci_latex(x)
    mantissa, exponent = split(@sprintf("%.3e", x), "e")
    exponent_int = parse(Int, exponent)
    return "\$$(mantissa)\\times10^{$(exponent_int)}\$"
end


end # module
