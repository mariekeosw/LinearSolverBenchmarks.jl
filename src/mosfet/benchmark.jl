# Benchmark Setup for MOSFET simulation: Runtimes Setup & Algorithms Setup

################################################################################
# Runtimes Setup
################################################################################
# Parameters of Interest: for one specific λ1, V_G or V_D
Base.@kwdef struct BenchmarkStep
    value::Float64 = 0.0     # λ1, V_G, V_D
    newton_iters::Int = 0

    tlinsolve::Float64 = 0.0
    tlinsolve_setup::Float64 = 0.0
    tlinsolve_solve::Float64 = 0.0
end

# Parameters of Interest: for all λ1, V_G or V_D values
# This struct is used for the detailed latex table (see below)
Base.@kwdef mutable struct BenchmarkSeries
    equilibrium::Vector{BenchmarkStep} = BenchmarkStep[]
    gate::Vector{BenchmarkStep} = BenchmarkStep[]
    drain::Vector{BenchmarkStep} = BenchmarkStep[]
end

# Runtimes of Interest: sum over all "λ1, V_G or V_D runtimes"
# This struct is used for the actual benchmark table (see benchmark() and PostProcess/Tab1Benchmark2DGrid.jl)
Base.@kwdef mutable struct BenchmarkSummary
    newton_iters::Int = 0

    tlinsolve::Float64 = 0.0
    tlinsolve_setup::Float64 = 0.0
    tlinsolve_solve::Float64 = 0.0
end

function sum_runtimes_mean_newton_iters(benchmark_data::BenchmarkSeries)

    newton_iters_equilibrium = sum(ni.newton_iters for ni in benchmark_data.equilibrium)
    newton_iters_equilibrium_length = length(ni.newton_iters for ni in benchmark_data.equilibrium)
    newton_iters_gate = sum(ni.newton_iters for ni in benchmark_data.gate)
    newton_iters_gate_length = length(ni.newton_iters for ni in benchmark_data.gate)
    newton_iters_drain = sum(ni.newton_iters for ni in benchmark_data.drain)
    newton_iters_drain_length = length(ni.newton_iters for ni in benchmark_data.drain)

    tlinsolve_equilibrium = sum(rt.tlinsolve for rt in benchmark_data.equilibrium)
    tlinsolve_gate = sum(rt.tlinsolve for rt in benchmark_data.gate)
    tlinsolve_drain = sum(rt.tlinsolve for rt in benchmark_data.drain)

    tlinsolve_setup_equilibrium = sum(rt.tlinsolve_setup for rt in benchmark_data.equilibrium)
    tlinsolve_setup_gate = sum(rt.tlinsolve_setup for rt in benchmark_data.gate)
    tlinsolve_setup_drain = sum(rt.tlinsolve_setup for rt in benchmark_data.drain)

    tlinsolve_solve_equilibrium = sum(rt.tlinsolve_solve for rt in benchmark_data.equilibrium)
    tlinsolve_solve_gate = sum(rt.tlinsolve_solve for rt in benchmark_data.gate)
    tlinsolve_solve_drain = sum(rt.tlinsolve_solve for rt in benchmark_data.drain)

    return BenchmarkSummary(
        newton_iters = round(Int, (newton_iters_equilibrium + newton_iters_gate + newton_iters_drain) / (newton_iters_equilibrium_length + newton_iters_gate_length + newton_iters_drain_length)),
        tlinsolve = tlinsolve_equilibrium + tlinsolve_gate + tlinsolve_drain,
        tlinsolve_setup = tlinsolve_setup_equilibrium + tlinsolve_setup_gate + tlinsolve_setup_drain,
        tlinsolve_solve = tlinsolve_solve_equilibrium + tlinsolve_solve_gate + tlinsolve_solve_drain
    )
end


################################################################################
# Algorithms Setup
################################################################################
Base.@kwdef struct SolverConfig
    label::String
    alg
end

################################################################################
# Comparison of Direct Solvers
################################################################################
# All direct solvers considered in the 2D analysis
function solver_configs_direct_all()
    return (
        SolverConfig("UMFPACK", UMFPACKFactorization()),            # LinearSolve.jl
        SolverConfig("KLU", KLUFactorization()),                    # LinearSolve.jl
        SolverConfig("MKLPardiso", MKLPardisoFactorize()),          # Pardiso.jl
        SolverConfig("Sparspak", SparspakFactorization()),          # LinearSolve.jl
    )
end

# MKLPardiso reference solver used for the 3D analysis
function solver_configs_direct_pardiso()
    return (
        SolverConfig("MKLPardiso", MKLPardisoFactorize()),          # Pardiso.jl
    )
end

# Reduced direct-solver configuration excluding KLU
function solver_configs_direct_without_klu()
    return (
        SolverConfig("UMFPACK", UMFPACKFactorization()),            # LinearSolve.jl
        SolverConfig("MKLPardiso", MKLPardisoFactorize()),          # Pardiso.jl
        SolverConfig("Sparspak", SparspakFactorization()),          # LinearSolve.jl
    )
end

################################################################################
# Comparison of Equation-Based and Node-Based Partitioning
################################################################################
function solver_configs_gmres_ilu_partitioning_comparison()
    # Equation-based block Jacobi (EBJ): one block for each physical unknown, containing its values at all nodes
    equation_partitioning = A -> [s:3:size(A, 1) for s in 1:3]

    ilut_droptol001_ebj = BlockPreconBuilder(
        precs = ExtendableSparse.ILUTPreconBuilder(droptol = 0.01),
        partitioning = equation_partitioning,
    )

    ilut_droptol01_ebj = BlockPreconBuilder(
        precs = ExtendableSparse.ILUTPreconBuilder(droptol = 0.1),
        partitioning = equation_partitioning,
    )

    ilu0_ebj = BlockPreconBuilder(
        precs = ExtendableSparse.ILUZeroPreconBuilder(),
        partitioning = equation_partitioning,
    )

    # Node-based block Jacobi (NBJ): one 3 × 3 diagonal block for the three unknowns at each grid node
    node_partitioning = A -> [((k - 1) * 3 + 1):(k * 3) for k in 1:(size(A, 1) ÷ 3)]

    ilut_droptol001_nbj = BlockPreconBuilder(
        precs = ExtendableSparse.ILUTPreconBuilder(droptol = 0.01),
        partitioning = node_partitioning,
    )

    ilut_droptol01_nbj = BlockPreconBuilder(
        precs = ExtendableSparse.ILUTPreconBuilder(droptol = 0.1),
        partitioning = node_partitioning,
    )

    ilu0_nbj = BlockPreconBuilder(
        precs = ExtendableSparse.ILUZeroPreconBuilder(),
        partitioning = node_partitioning,
    )

    # Block ILU(0) with blocksize = 3 was excluded because it did not converge
    # for the equilibrium problem within the prescribed iteration limit.

    return (
        SolverConfig("MKLPardiso", MKLPardisoFactorize()),
        SolverConfig("GMRES-EBJ-ILUT-0.01", KrylovJL_GMRES(precs = ilut_droptol001_ebj)),
        SolverConfig("GMRES-EBJ-ILUT-0.1", KrylovJL_GMRES(precs = ilut_droptol01_ebj)),
        SolverConfig("GMRES-EBJ-ILU0", KrylovJL_GMRES(precs = ilu0_ebj)),
        SolverConfig("GMRES-NBJ-ILUT-0.01", KrylovJL_GMRES(precs = ilut_droptol001_nbj)),
        SolverConfig("GMRES-NBJ-ILUT-0.1", KrylovJL_GMRES(precs = ilut_droptol01_nbj)),
        SolverConfig("GMRES-NBJ-ILU0", KrylovJL_GMRES(precs = ilu0_nbj)),
    )
end

function solver_configs_gmres_amg_partitioning_comparison()
    # Equation-based block Jacobi (EBJ): one block for each physical unknown, containing its values at all nodes
    equation_partitioning = A -> [s:3:size(A, 1) for s in 1:3]

    amgcl_sa_spai0_ebj = BlockPreconBuilder(
        precs = AMGCLWrap.AMGPreconBuilder(coarsening = AMGCLWrap.SmoothedAggregationCoarsening(), relax = AMGCLWrap.SPAI0Relaxation()),
        partitioning = equation_partitioning,
    )

    amg_rs_ebj = BlockPreconBuilder(
        precs = AlgebraicMultigrid.RugeStubenPreconBuilder(),
        partitioning = equation_partitioning,
    )

    # Node-based block Jacobi (NBJ): one 3 × 3 diagonal block for the three unknowns at each grid node
    node_partitioning = A -> [((k - 1) * 3 + 1):(k * 3) for k in 1:(size(A, 1) ÷ 3)]

    amgcl_sa_spai0_nbj = BlockPreconBuilder(
        precs = AMGCLWrap.AMGPreconBuilder(coarsening = AMGCLWrap.SmoothedAggregationCoarsening(), relax = AMGCLWrap.SPAI0Relaxation()),
        partitioning = node_partitioning,
    )

    amg_rs_nbj = BlockPreconBuilder(
        precs = AlgebraicMultigrid.RugeStubenPreconBuilder(),
        partitioning = node_partitioning,
    )

    # Block-AMG variants with blocksize = 3 were excluded:
    # AMGCL produced inconsistent results or convergence failures, while
    # AlgebraicMultigrid raised a linear solver error.

    return (
        SolverConfig("MKLPardiso", MKLPardisoFactorize()),
        SolverConfig("GMRES-EBJ-AMGCL-SA-SPAI0", KrylovJL_GMRES(precs = amgcl_sa_spai0_ebj)),
        SolverConfig("GMRES-EBJ-AMG-RS", KrylovJL_GMRES(precs = amg_rs_ebj)),
        SolverConfig("GMRES-NBJ-AMGCL-SA-SPAI0", KrylovJL_GMRES(precs = amgcl_sa_spai0_nbj)),
        SolverConfig("GMRES-NBJ-AMG-RS", KrylovJL_GMRES(precs = amg_rs_nbj)),
    )
end

function solver_configs_bicgstab_ilu_partitioning_comparison()
    # Equation-based block Jacobi (EBJ): one block for each physical unknown, containing its values at all nodes
    equation_partitioning = A -> [s:3:size(A, 1) for s in 1:3]

    ilut_droptol001_ebj = BlockPreconBuilder(
        precs = ExtendableSparse.ILUTPreconBuilder(droptol = 0.01),
        partitioning = equation_partitioning,
    )

    ilut_droptol01_ebj = BlockPreconBuilder(
        precs = ExtendableSparse.ILUTPreconBuilder(droptol = 0.1),
        partitioning = equation_partitioning,
    )

    ilu0_ebj = BlockPreconBuilder(
        precs = ExtendableSparse.ILUZeroPreconBuilder(),
        partitioning = equation_partitioning,
    )

    # Node-based block Jacobi (NBJ): one 3 × 3 diagonal block for the three unknowns at each grid node
    node_partitioning = A -> [((k - 1) * 3 + 1):(k * 3) for k in 1:(size(A, 1) ÷ 3)]

    ilut_droptol001_nbj = BlockPreconBuilder(
        precs = ExtendableSparse.ILUTPreconBuilder(droptol = 0.01),
        partitioning = node_partitioning,
    )

    ilut_droptol01_nbj = BlockPreconBuilder(
        precs = ExtendableSparse.ILUTPreconBuilder(droptol = 0.1),
        partitioning = node_partitioning,
    )

    ilu0_nbj = BlockPreconBuilder(
        precs = ExtendableSparse.ILUZeroPreconBuilder(),
        partitioning = node_partitioning,
    )

    # Block ILU(0) with blocksize = 3 was excluded because it did not converge
    # for the equilibrium problem within the prescribed iteration limit.

    return (
        SolverConfig("MKLPardiso", MKLPardisoFactorize()),
        SolverConfig("BiCGSTAB-EBJ-ILUT-0.01", KrylovJL_BICGSTAB(precs = ilut_droptol001_ebj)),
        SolverConfig("BiCGSTAB-EBJ-ILUT-0.1", KrylovJL_BICGSTAB(precs = ilut_droptol01_ebj)),
        SolverConfig("BiCGSTAB-EBJ-ILU0", KrylovJL_BICGSTAB(precs = ilu0_ebj)),
        SolverConfig("BiCGSTAB-NBJ-ILUT-0.01", KrylovJL_BICGSTAB(precs = ilut_droptol001_nbj)),
        SolverConfig("BiCGSTAB-NBJ-ILUT-0.1", KrylovJL_BICGSTAB(precs = ilut_droptol01_nbj)),
        SolverConfig("BiCGSTAB-NBJ-ILU0", KrylovJL_BICGSTAB(precs = ilu0_nbj)),
    )
end

function solver_configs_bicgstab_amg_partitioning_comparison()
    # Equation-based block Jacobi (EBJ): one block for each physical unknown, containing its values at all nodes
    equation_partitioning = A -> [s:3:size(A, 1) for s in 1:3]

    amgcl_sa_spai0_ebj = BlockPreconBuilder(
        precs = AMGCLWrap.AMGPreconBuilder(coarsening = AMGCLWrap.SmoothedAggregationCoarsening(), relax = AMGCLWrap.SPAI0Relaxation()),
        partitioning = equation_partitioning,
    )

    # Node-based block Jacobi (NBJ): one 3 × 3 diagonal block for the three unknowns at each grid node
    node_partitioning = A -> [((k - 1) * 3 + 1):(k * 3) for k in 1:(size(A, 1) ÷ 3)]

    amgcl_sa_spai0_nbj = BlockPreconBuilder(
        precs = AMGCLWrap.AMGPreconBuilder(coarsening = AMGCLWrap.SmoothedAggregationCoarsening(), relax = AMGCLWrap.SPAI0Relaxation()),
        partitioning = node_partitioning,
    )

    # Block-AMG variants with blocksize = 3 were excluded because they produced
    # inconsistent results or solver failures for the considered MOSFET problem.

    return (
        SolverConfig("MKLPardiso", MKLPardisoFactorize()),
        SolverConfig("BiCGSTAB-EBJ-AMGCL-SA-SPAI0", KrylovJL_BICGSTAB(precs = amgcl_sa_spai0_ebj)),
        SolverConfig("BiCGSTAB-NBJ-AMGCL-SA-SPAI0", KrylovJL_BICGSTAB(precs = amgcl_sa_spai0_nbj)),
    )
end

################################################################################
# Compare different Preconditioners
################################################################################
# GMRES with equation-based block-Jacobi (EBJ) preconditioning
function solver_configs_gmres_ebj_preconditioner_comparison()
    equation_partitioning = A -> [s:3:size(A, 1) for s in 1:3]

    ilut_droptol001 = BlockPreconBuilder(
        precs = ExtendableSparse.ILUTPreconBuilder(droptol = 0.01),
        partitioning = equation_partitioning,
    )

    ilut_droptol01 = BlockPreconBuilder(
        precs = ExtendableSparse.ILUTPreconBuilder(droptol = 0.1),
        partitioning = equation_partitioning,
    )

    ilu0 = BlockPreconBuilder(
        precs = ExtendableSparse.ILUZeroPreconBuilder(),
        partitioning = equation_partitioning,
    )

    amgcl_sa_spai0 = BlockPreconBuilder(
        precs = AMGCLWrap.AMGPreconBuilder(coarsening = AMGCLWrap.SmoothedAggregationCoarsening(), relax = AMGCLWrap.SPAI0Relaxation()),
        partitioning = equation_partitioning,
    )

    amgcl_rs_spai0 = BlockPreconBuilder(
        precs = AMGCLWrap.AMGPreconBuilder(coarsening = AMGCLWrap.RugeStubenCoarsening(), relax = AMGCLWrap.SPAI0Relaxation()),
        partitioning = equation_partitioning,
    )

    amgcl_sa_ilu0 = BlockPreconBuilder(
        precs = AMGCLWrap.AMGPreconBuilder(coarsening = AMGCLWrap.SmoothedAggregationCoarsening(), relax = AMGCLWrap.ILU0Relaxation()),
        partitioning = equation_partitioning,
    )

    amgcl_rs_ilu0 = BlockPreconBuilder(
        precs = AMGCLWrap.AMGPreconBuilder(coarsening = AMGCLWrap.RugeStubenCoarsening(), relax = AMGCLWrap.ILU0Relaxation()),
        partitioning = equation_partitioning,
    )

    # AlgebraicMultigrid.jl variants were excluded because they either failed
    # solution validation or raised solver errors for the considered MOSFET problem.

    return (
        SolverConfig("MKLPardiso", MKLPardisoFactorize()),
        SolverConfig("GMRES-EBJ-ILUT-0.01", KrylovJL_GMRES(precs = ilut_droptol001)),
        SolverConfig("GMRES-EBJ-ILUT-0.1", KrylovJL_GMRES(precs = ilut_droptol01)),
        SolverConfig("GMRES-EBJ-ILU0", KrylovJL_GMRES(precs = ilu0)),
        SolverConfig("GMRES-EBJ-AMGCL-SA-SPAI0", KrylovJL_GMRES(precs = amgcl_sa_spai0)),
        SolverConfig("GMRES-EBJ-AMGCL-RS-SPAI0", KrylovJL_GMRES(precs = amgcl_rs_spai0)),
        SolverConfig("GMRES-EBJ-AMGCL-SA-ILU0", KrylovJL_GMRES(precs = amgcl_sa_ilu0)),
        SolverConfig("GMRES-EBJ-AMGCL-RS-ILU0", KrylovJL_GMRES(precs = amgcl_rs_ilu0)),
    )
end

# BiCGSTAB with equation-based block-Jacobi (EBJ) preconditioning
function solver_configs_bicgstab_ebj_preconditioner_comparison()
    equation_partitioning = A -> [s:3:size(A, 1) for s in 1:3]

    ilut_droptol001 = BlockPreconBuilder(
        precs = ExtendableSparse.ILUTPreconBuilder(droptol = 0.01),
        partitioning = equation_partitioning,
    )

    ilut_droptol01 = BlockPreconBuilder(
        precs = ExtendableSparse.ILUTPreconBuilder(droptol = 0.1),
        partitioning = equation_partitioning,
    )

    ilu0 = BlockPreconBuilder(
        precs = ExtendableSparse.ILUZeroPreconBuilder(),
        partitioning = equation_partitioning,
    )

    amgcl_sa_spai0 = BlockPreconBuilder(
        precs = AMGCLWrap.AMGPreconBuilder(coarsening = AMGCLWrap.SmoothedAggregationCoarsening(), relax = AMGCLWrap.SPAI0Relaxation()),
        partitioning = equation_partitioning,
    )

    amgcl_rs_spai0 = BlockPreconBuilder(
        precs = AMGCLWrap.AMGPreconBuilder(coarsening = AMGCLWrap.RugeStubenCoarsening(), relax = AMGCLWrap.SPAI0Relaxation()),
        partitioning = equation_partitioning,
    )

    amgcl_sa_ilu0 = BlockPreconBuilder(
        precs = AMGCLWrap.AMGPreconBuilder(coarsening = AMGCLWrap.SmoothedAggregationCoarsening(), relax = AMGCLWrap.ILU0Relaxation()),
        partitioning = equation_partitioning,
    )

    amgcl_rs_ilu0 = BlockPreconBuilder(
        precs = AMGCLWrap.AMGPreconBuilder(coarsening = AMGCLWrap.RugeStubenCoarsening(), relax = AMGCLWrap.ILU0Relaxation()),
        partitioning = equation_partitioning,
    )

    # AlgebraicMultigrid.jl variants were excluded because they either failed
    # solution validation or raised solver errors for the considered MOSFET problem.

    return (
        SolverConfig("MKLPardiso", MKLPardisoFactorize()),
        SolverConfig("BiCGSTAB-EBJ-ILUT-0.01", KrylovJL_BICGSTAB(precs = ilut_droptol001)),
        SolverConfig("BiCGSTAB-EBJ-ILUT-0.1", KrylovJL_BICGSTAB(precs = ilut_droptol01)),
        SolverConfig("BiCGSTAB-EBJ-ILU0", KrylovJL_BICGSTAB(precs = ilu0)),
        SolverConfig("BiCGSTAB-EBJ-AMGCL-SA-SPAI0", KrylovJL_BICGSTAB(precs = amgcl_sa_spai0)),
        SolverConfig("BiCGSTAB-EBJ-AMGCL-RS-SPAI0", KrylovJL_BICGSTAB(precs = amgcl_rs_spai0)),
        SolverConfig("BiCGSTAB-EBJ-AMGCL-SA-ILU0", KrylovJL_BICGSTAB(precs = amgcl_sa_ilu0)),
        SolverConfig("BiCGSTAB-EBJ-AMGCL-RS-ILU0", KrylovJL_BICGSTAB(precs = amgcl_rs_ilu0)),
    )
end

################################################################################
# Selected Solver Configurations for Larger Problem Sizes
################################################################################
function solver_configs_gmres_bicgstab_ebj()
    equation_partitioning = A -> [s:3:size(A, 1) for s in 1:3]

    ilu0_ebj = BlockPreconBuilder(
        precs = ExtendableSparse.ILUZeroPreconBuilder(),
        partitioning = equation_partitioning,
    )

    amgcl_sa_spai0_ebj = BlockPreconBuilder(
        precs = AMGCLWrap.AMGPreconBuilder(coarsening = AMGCLWrap.SmoothedAggregationCoarsening(), relax = AMGCLWrap.SPAI0Relaxation()),
        partitioning = equation_partitioning,
    )

    return (
        SolverConfig("MKLPardiso", MKLPardisoFactorize()),
        SolverConfig("GMRES-EBJ-ILU0", KrylovJL_GMRES(precs = ilu0_ebj)),
        SolverConfig("GMRES-EBJ-AMGCL-SA-SPAI0", KrylovJL_GMRES(precs = amgcl_sa_spai0_ebj)),
        SolverConfig("BiCGSTAB-EBJ-ILU0", KrylovJL_BICGSTAB(precs = ilu0_ebj)),
        SolverConfig("BiCGSTAB-EBJ-AMGCL-SA-SPAI0", KrylovJL_BICGSTAB(precs = amgcl_sa_spai0_ebj)),
    )
end

function solver_configs_gmres_bicgstab_ebj_without_direct()
    equation_partitioning = A -> [s:3:size(A, 1) for s in 1:3]

    ilu0_ebj = BlockPreconBuilder(
        precs = ExtendableSparse.ILUZeroPreconBuilder(),
        partitioning = equation_partitioning,
    )

    amgcl_sa_spai0_ebj = BlockPreconBuilder(
        precs = AMGCLWrap.AMGPreconBuilder(coarsening = AMGCLWrap.SmoothedAggregationCoarsening(), relax = AMGCLWrap.SPAI0Relaxation()),
        partitioning = equation_partitioning,
    )

    return (
        SolverConfig("GMRES-EBJ-ILU0", KrylovJL_GMRES(precs = ilu0_ebj)),
        SolverConfig("GMRES-EBJ-AMGCL-SA-SPAI0", KrylovJL_GMRES(precs = amgcl_sa_spai0_ebj)),
        SolverConfig("BiCGSTAB-EBJ-ILU0", KrylovJL_BICGSTAB(precs = ilu0_ebj)),
        SolverConfig("BiCGSTAB-EBJ-AMGCL-SA-SPAI0", KrylovJL_BICGSTAB(precs = amgcl_sa_spai0_ebj)),
    )
end
