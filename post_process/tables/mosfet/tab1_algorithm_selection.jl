module Tab1AlgorithmSelection

using JLD2
using Printf

# Helper function: find test entry
function get_values(tests, tests_phy, methodname; ref_id = 1)

    # Falls der Schlüssel nicht exisitiert - leerer Eintrag
    if !haskey(tests, methodname) || !haskey(tests_phy, methodname)
        return ("", "")
    end

    test = tests[methodname][ref_id]
    test_phy = tests_phy[methodname][ref_id]

    return (
        @sprintf("%.3e", test),
        @sprintf("%.3e", test_phy),
    )
end

# How to use this function?
# 1) includet("post_process/tables/mosfet/tab1_algorithm_selection.jl")
# 2) Tab1AlgorithmSelection.export_benchmark_to_latex("xxx.jld2", "xxx.tex")
function export_benchmark_to_latex(filename_in, filename_out)

    @load filename_in labels tests tests_phy

    # Erste Tabellenzeile
    methods1 = ["MKLPardiso", "BiLQ-ILU0", "QMR-ILU0", "CGS-ILU0", "BiCGSTAB-ILU0"]

    # Zweite Tabellenzeile
    methods2 = ["DIOM-ILU0", "FOM-ILU0", "DQGMRES-ILU0", "GMRES-ILU0", "FGMRES-ILU0"]

    open(filename_out, "w") do io

        println(io, "\\begin{table}[h]")
        println(io, "\\centering")
        println(io, "\\footnotesize")
        println(io, "\\caption{Benchmark results for Krylov methods using Block Jacobi preconditioning with ILU(0)}")
        println(io, "\\begin{tabular}{lccccc}")
        println(io, "\\toprule")

        #######################################################################
        # First block
        #######################################################################

        println(io, "Method & ", join(methods1, " & "), " \\\\")

        println(io, "\\midrule")

        testvals1 = [get_values(tests, tests_phy, m)[1] for m in methods1]
        testphy1 = [get_values(tests, tests_phy, m)[2] for m in methods1]

        println(io, "Test 1 & ", join(testvals1, " & "), " \\\\[0.1cm]")

        println(io, "Test 2 & ", join(testphy1, " & "), " \\\\")

        println(io, "\\bottomrule")
        println(io, "\\addlinespace[1.0em]")
        println(io, "\\toprule")

        #######################################################################
        # Second block
        #######################################################################

        println(io, "Method & ", join(methods2, " & "), " \\\\")

        println(io, "\\midrule")

        testvals2 = [get_values(tests, tests_phy, m)[1] for m in methods2]
        testphy2 = [get_values(tests, tests_phy, m)[2] for m in methods2]

        println(io, "Test 1 & ", join(testvals2, " & "), " \\\\[0.1cm]")

        println(io, "Test 2 & ", join(testphy2, " & "), " \\\\")

        println(io, "\\bottomrule")
        println(io, "\\end{tabular}")
        println(io, "\\label{tab:overview_krylov}")
        println(io, "\\end{table}")
    end

    return nothing
end # function

end # module
