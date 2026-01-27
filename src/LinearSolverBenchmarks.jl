module LinearSolverBenchmarks

module ModelProblem
    using LinearAlgebra: I, lu
    using LinearSolve
    using Pardiso
    using Sparspak
    using ILUZero: ilu0
    using IncompleteLU
    using AMGCLWrap
    using AlgebraicMultigrid
    using LinearSolvePyAMG
    using Krylov

    include("model_problem/benchmark.jl")
    export SolverConfig
    export solver_configs_direct
    export solver_configs_direct_pardiso
    export sym_solver_configs_cg
    export sym_solver_configs_cg_ilu
    export sym_solver_configs_cg_amgcl
    export sym_solver_configs_cg_amg
    export sym_solver_configs_cg_pyamg
    export nsym_solver_configs_krylov
    export nsym_solver_configs_cgs_bicgstab
    export nsym_solver_configs_gmres
    export nsym_solver_configs_krylov_ilu
    export nsym_solver_configs_krylov_amgcl_spai0
    export nsym_solver_configs_krylov_amgcl_ilu0
    export nsym_solver_configs_krylov_amg
    export nsym_solver_configs_krylov_pyamg

end # module ModelProblem

module Mosfet
    using ChargeTransport
    using ExtendableGrids
    using VoronoiFVM
    using ExtendableSparse
    using Krylov
    using LinearSolve
    using LinearSolvePyAMG
    using AMGCLWrap
    using AlgebraicMultigrid
    using Pardiso
    using LaTeXStrings
    using Printf

    include("mosfet/benchmark.jl")
    export BenchmarkStep
    export BenchmarkSeries
    export BenchmarkSummary
    export sum_runtimes_mean_newton_iters
    export SolverConfig
    export solver_configs_direct_all
    export solver_configs_direct_pardiso
    export solver_configs_direct_without_klu
    export solver_configs_gmres_ilu_partitioning_comparison
    export solver_configs_gmres_amg_partitioning_comparison
    export solver_configs_bicgstab_ilu_partitioning_comparison
    export solver_configs_bicgstab_amg_partitioning_comparison
    export solver_configs_gmres_ebj_preconditioner_comparison
    export solver_configs_bicgstab_ebj_preconditioner_comparison
    export solver_configs_gmres_bicgstab_ebj
    export solver_configs_gmres_bicgstab_ebj_without_direct

    include("mosfet/export.jl")
    export tridata
    export plot_heatmap_matrix
    export export_detailed_benchmark_to_latex
    export plot_region_doping_2D
    export plot_nodal_doping_2D
    export plotsurface_nodal_doping_2D


    include("mosfet/grid.jl")
    export generate_2Dgrid
    export generate_3Dgrid
    export generate_3Dgrid_half

    include("mosfet/solve.jl")
    export equilibrium_solve!
    export gate_solve!
    export drain_solve!

    parametersdir(args...) = joinpath(@__DIR__, "mosfet", "setups", args...)

    for parameter_set in [
            :tesca_setup,
        ]
        include(parametersdir("$(parameter_set).jl"))
        @eval export $parameter_set
    end
end # module Mosfet

end # module LinearSolverBenchmarks
