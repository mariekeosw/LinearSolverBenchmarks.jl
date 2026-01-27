#=
MOSFET Simulation and Linear Solver Benchmark

This script performs stationary drift-diffusion simulations for an n-channel
MOSFET consisting of a single material, namely silicon.

Two doping concepts are used in the overall framework. This script
uses region-wise doping. For the alternative node-wise doping approach, see
"script/analyse_mosfet_doping.jl".

The file is divided into two main parts:
1. the actual semiconductor simulation, implemented in main(), and
2. the benchmarking routine benchmark(), which compares different linear solver
configurations for a sequence of grid refinements.
=#

module MosfetBenchmarkRun

using ChargeTransport
using VoronoiFVM          # for IV Curve, Current Density Plot
using ExtendableGrids
using LinearSolverBenchmarks: Mosfet
using JLD2

################################################################################
# Semiconductor Simulation
################################################################################
function main(;
        Documentation = true,
        ###############################
        grid_dimension = 2,                 # 2D or 3D
        grid_refinement = 0,                # 2D: 0 with 930 nodes, 1 with 3599 nodes, ... and 3D: 0 with 28830 nodes, 1 with 219539 nodes, ...
        ###############################
        initial_value = :charge_neutrality, # initial value for the equilibrium run (either: :charge_neutrality or :lambda_embedding)
        strategy = nothing,                 # Default for sparse matrices is UMFPACKFactorization() (see VoronoiFVM.jl/src/vfvm_solver.jl)
        ###############################
        parameter_set = Mosfet.tesca_setup,
        Plotter = nothing                   # Heatmap
    )

    # Initialize storage for benchmark data (runtimes, iterations)
    benchmark_data = Mosfet.BenchmarkSeries()

    ################################################################################
    if Documentation
        println("Set up grid and regions")
    end
    ################################################################################

    # Select grid based on grid_dimension
    if grid_dimension == 2
        grid = Mosfet.generate_2Dgrid(grid_refinement)
    elseif grid_dimension == 3
        grid = Mosfet.generate_3Dgrid_half(grid_refinement)
    else
        error("Please select grid dimension 2 or 3.")
    end

    # Store number of grid nodes (from that the Dimension of the Jacobian can be derived)
    num_nodes = size(grid[Coordinates], 2)

    if Documentation
        println("*** done\n")
    end

    ################################################################################
    if Documentation
        println("Define physical parameters and model")
    end
    ################################################################################

    @local_unitfactors μm cm s ns V K

    constants = ChargeTransport.teSCA_constants

    # Import parameters from parameter file (TeSCA parameter set)
    p = parameter_set()

    # Initialize Data instance
    data = Data(grid, p.numberOfCarriers, constants = constants)

    data.modelType = Stationary

    data.F .= Boltzmann

    data.bulkRecombination = set_bulk_recombination(;
        iphin = p.iphin, iphip = p.iphip,
        bulk_recomb_Auger = true,
        bulk_recomb_radiative = false,
        bulk_recomb_SRH = true,
    )

    data.boundaryType[p.bregion_gate] = GateContact
    data.boundaryType[p.bregion_drain] = OhmicContact
    data.boundaryType[p.bregion_source] = OhmicContact
    data.boundaryType[p.bregion_bulk] = OhmicContact

    # Flux discretization - depends on statistic (Boltzmann -> Scharfetter Gummel)
    data.fluxApproximation .= ScharfetterGummel

    data.params = ChargeTransport.Params(p)

    # Definition ChargeTransport System
    ctsys = ChargeTransport.System(grid, data, unknown_storage = :sparse)

    if Documentation
        ChargeTransport.show_params(ctsys)
        println("*** done\n")
    end

    ################################################################################
    if Documentation
        println("Define control parameters for Solver")
    end
    ################################################################################

    control = ChargeTransport.SolverControl()
    control.verbose = false
    control.log = true
    # Newton control parameters
    control.maxiters = 20
    control.abstol = 1.0e-10
    control.reltol = 1.0e-10
    control.tol_round = 1.0e-10
    control.damp_initial = 0.5
    control.max_round = 3
    # Iterative linear control parameters
    control.method_linear = strategy
    control.reltol_linear = 1.0e-4
    control.abstol_linear = 1.0e-8
    control.maxiters_linear = 100
    # Embedding
    control.Δp = 1.0e-20
    control.Δp_min = 1.0e-20
    control.Δp_max = 1.0
    control.Δp_grow = 10
    control.Δp_decrease = 0.5
    control.Δu_opt = 2.0
    control.handle_exceptions = true

    if Documentation
        println("*** done\n")
    end

    ################################################################################
    if Documentation
        println("Compute solution in thermodynamic equilibrium")
    end
    ################################################################################

    # Solution in Equilibrium
    # Argument "benchmark_data" is passed to equilibrium_solve! to store runtimes, number of iterations
    solution_eq = Mosfet.equilibrium_solve!(ctsys, initial_value, grid_dimension, control = control, benchmark_data = benchmark_data)

    # Plot Heatmap Equilibrium
    #heatmap_eq = Mosfet.plot_heatmap_matrix(ctsys.fvmsys, solution_eq; blocksize = 30, permutation = false, Plotter = Plotter)

    if Documentation
        println("*** done\n")
    end

    ################################################################################
    if Documentation
        println("Bias loop")
    end
    ################################################################################

    # Adjust control-embedding parameters for full solution
    control.Δp = 0.2
    control.Δp_min = 1.0e-3
    control.Δp_max = 1.0
    control.Δp_grow = 1.2
    control.Δp_decrease = 0.5
    control.Δu_opt = 1.5
    control.handle_exceptions = true

    inival = copy(solution_eq)
    solution = copy(solution_eq)

    # Gate bias: start at 0 V and increase up to 5 V
    biasValue_gate_range = [0.0, 5.0]
    solution = Mosfet.gate_solve!(ctsys, biasValue_gate_range, inival = inival, control = control, p = p, benchmark_data = benchmark_data)
    inival .= solution

    # Drain bias: start at 0 V and increase up to 5 V
    biasValue_drain_range = [0.0, 5.0]
    solution, biasValues_drain, IV_drain = Mosfet.drain_solve!(ctsys, biasValue_drain_range, grid_dimension, inival = inival, control = control, p = p, benchmark_data = benchmark_data)

    # Plot Heatmap Full Solution
    #heatmap = Mosfet.plot_heatmap_matrix(ctsys.fvmsys, solution; blocksize = 30, permutation = false, Plotter = Plotter)

    if Documentation
        println("*** done\n")
    end

    # BenchmarkSummary
    # (1) mean value of newton iterations
    # (2) total runtimes for equilibrium, gate bias and drain bias (sum)
    benchmark_summary = Mosfet.sum_runtimes_mean_newton_iters(benchmark_data)

    # Unphysical Testvalue
    testval = sum(solution) / length(solution)

    # Physical Testvalue - IV value for drainvoltage = 5V
    testval_phy = IV_drain[end]

    # Export detailed iteration statistics (Newton steps, i.e. number of linear systems solved)
    #Mosfet.export_detailed_benchmark_to_latex(benchmark_data; filename_out = "linear_systems_counter_le.tex")

    # Save simulation results for post-processing and visualization
    @save "data_mosfet_$(grid_dimension)d_ref$(grid_refinement).jld2" grid_refinement solution_eq solution biasValues_drain IV_drain

    return num_nodes, testval, testval_phy, benchmark_summary.newton_iters, benchmark_summary.tlinsolve, benchmark_summary.tlinsolve_setup, benchmark_summary.tlinsolve_solve
end # function main


################################################################################
# Benchmark the Semiconductor Simulation
################################################################################
# Call function: benchmark("test.jld2")
function benchmark(
        filename;
        grid_dimension = 2,
        grid_refinement_max = 0,
        initial_value = :charge_neutrality,
        solver_config = :direct_all
    )

    if solver_config == :direct_all
        configs = Mosfet.solver_configs_direct_all()

    elseif solver_config == :direct_pardiso
        configs = Mosfet.solver_configs_direct_pardiso()

    elseif solver_config == :direct_all_without_klu
        configs = Mosfet.solver_configs_direct_without_klu()

    elseif solver_config == :gmres_ilu_part_comparison
        configs = Mosfet.solver_configs_gmres_ilu_partitioning_comparison()

    elseif solver_config == :gmres_amg_part_comparison
        configs = Mosfet.solver_configs_gmres_amg_partitioning_comparison()

    elseif solver_config == :bicgstab_ilu_part_comparison
        configs = Mosfet.solver_configs_bicgstab_ilu_partitioning_comparison()

    elseif solver_config == :bicgstab_amg_part_comparison
        configs = Mosfet.solver_configs_bicgstab_amg_partitioning_comparison()

    elseif solver_config == :gmres_ebj
        configs = Mosfet.solver_configs_gmres_ebj_preconditioner_comparison()

    elseif solver_config == :bicgstab_ebj
        configs = Mosfet.solver_configs_bicgstab_ebj_preconditioner_comparison()

    elseif solver_config == :gmres_bicgstab_ebj
        configs = Mosfet.solver_configs_gmres_bicgstab_ebj()

    elseif solver_config == :gmres_bicgstab_ebj_without_direct
        configs = Mosfet.solver_configs_gmres_bicgstab_ebj_without_direct()

    else
        error("Unknown solver configuration = $solver_config")
    end

    labels = [config.label for config in configs]

    # Different numbers of nodes, size of the corresponding jac matrix differs
    ref = collect(0:grid_refinement_max)
    #ref = grid_refinement_max

    # Create a storage for results
    number_nodes = Dict{String, Vector{Int}}()
    walltimes = Dict{String, Vector{Float64}}()
    storage_tlinsolve = Dict{String, Vector{Float64}}()
    storage_tlinsolve_setup = Dict{String, Vector{Float64}}()
    storage_tlinsolve_solve = Dict{String, Vector{Float64}}()
    storage_newton_iters = Dict{String, Vector{Int}}()
    tests = Dict{String, Vector{Float64}}()
    tests_phy = Dict{String, Vector{Float64}}()

    for config in configs
        label = config.label
        alg = config.alg

        println("algorithm = $label")

        number_nodes_values = Vector{Int}(undef, length(ref))
        walltime_values = Vector{Float64}(undef, length(ref))
        tlinsolve_values = Vector{Float64}(undef, length(ref))
        tlinsolve_setup_values = Vector{Float64}(undef, length(ref))
        tlinsolve_solve_values = Vector{Float64}(undef, length(ref))
        newton_iters_values = Vector{Int}(undef, length(ref))
        test_values = Vector{Float64}(undef, length(ref))
        test_values_phy = Vector{Float64}(undef, length(ref))


        for j in eachindex(ref)
            println("refinementfactor = $(ref[j])")

            bt = @elapsed begin
                num_nodes, testval, testval_phy, newton_iters, tlinsolve, tlinsolve_setup, tlinsolve_solve = main(
                    Documentation = false,
                    grid_dimension = grid_dimension,
                    grid_refinement = ref[j],
                    initial_value = initial_value,
                    strategy = alg
                )
            end

            number_nodes_values[j] = num_nodes
            walltime_values[j] = bt
            tlinsolve_values[j] = tlinsolve
            tlinsolve_setup_values[j] = tlinsolve_setup
            tlinsolve_solve_values[j] = tlinsolve_solve
            newton_iters_values[j] = newton_iters
            test_values[j] = testval
            test_values_phy[j] = testval_phy

            println("num_nodes = $(num_nodes)")
            println("walltime = $bt")
            println("tlinsolve = $tlinsolve")
            println("tlinsolve_setup = $tlinsolve_setup")
            println("tlinsolve_solve = $tlinsolve_solve")
            println("newton_iters = $(newton_iters)")
            println("testval = $(testval)")
            println("testval_phy = $(testval_phy)")
            println("*** done\n")
        end

        number_nodes[label] = number_nodes_values
        walltimes[label] = walltime_values
        storage_tlinsolve[label] = tlinsolve_values
        storage_tlinsolve_setup[label] = tlinsolve_setup_values
        storage_tlinsolve_solve[label] = tlinsolve_solve_values
        storage_newton_iters[label] = newton_iters_values
        tests[label] = test_values
        tests_phy[label] = test_values_phy
    end

    @save filename labels ref number_nodes walltimes storage_tlinsolve storage_tlinsolve_setup storage_tlinsolve_solve storage_newton_iters tests tests_phy

    return nothing
end # function benchmark


################################################################################
# Benchmark Setup Server
################################################################################
function several_benchmarks(; initial_value = :charge_neutrality)

    function run_benchmark(
            filename;
            grid_dimension,
            grid_refinement_max,
            solver_config
        )

        # Warm-up
        MosfetBenchmarkRun.benchmark("warmup.jld2"; grid_dimension = grid_dimension, grid_refinement_max = 0, initial_value = initial_value, solver_config = solver_config)

        # Actual benchmark
        return @elapsed MosfetBenchmarkRun.benchmark(filename; grid_dimension = grid_dimension, grid_refinement_max = grid_refinement_max, initial_value = initial_value, solver_config = solver_config)
    end

    ############################################################################
    # Two-dimensional MOSFET problem
    ############################################################################

    tdirect_all_2d = run_benchmark("mosfet_2d_direct_all_7.jld2"; grid_dimension = 2, grid_refinement_max = 7, solver_config = :direct_all)

    tgmres_ilu_part_2d = run_benchmark("mosfet_2d_gmres_ilu_part_comparison_5.jld2"; grid_dimension = 2, grid_refinement_max = 5, solver_config = :gmres_ilu_part_comparison)

    tgmres_amg_part_2d = run_benchmark("mosfet_2d_gmres_amg_part_comparison_0.jld2"; grid_dimension = 2, grid_refinement_max = 0, solver_config = :gmres_amg_part_comparison)

    tbicgstab_ilu_part_2d = run_benchmark("mosfet_2d_bicgstab_ilu_part_comparison_5.jld2"; grid_dimension = 2, grid_refinement_max = 5, solver_config = :bicgstab_ilu_part_comparison)

    tbicgstab_amg_part_2d = run_benchmark("mosfet_2d_bicgstab_amg_part_comparison_0.jld2"; grid_dimension = 2, grid_refinement_max = 0, solver_config = :bicgstab_amg_part_comparison)

    tgmres_ebj_2d = run_benchmark("mosfet_2d_gmres_ebj_2.jld2"; grid_dimension = 2, grid_refinement_max = 2, solver_config = :gmres_ebj)

    tbicgstab_ebj_2d = run_benchmark("mosfet_2d_bicgstab_ebj_2.jld2"; grid_dimension = 2, grid_refinement_max = 2, solver_config = :bicgstab_ebj)

    tgmres_bicgstab_ebj_2d = run_benchmark("mosfet_2d_gmres_bicgstab_ebj_19.jld2"; grid_dimension = 2, grid_refinement_max = 19, solver_config = :gmres_bicgstab_ebj)

    ############################################################################
    # Three-dimensional MOSFET problem
    ############################################################################

    tdirect_pardiso_3d = run_benchmark("mosfet_3d_direct_pardiso_0.jld2"; grid_dimension = 3, grid_refinement_max = 0, solver_config = :direct_pardiso)

    tgmres_bicgstab_ebj_3d = run_benchmark("mosfet_3d_gmres_bicgstab_ebj_1.jld2"; grid_dimension = 3, grid_refinement_max = 1, solver_config = :gmres_bicgstab_ebj)

    tgmres_bicgstab_ebj_without_direct_3d = run_benchmark("mosfet_3d_gmres_bicgstab_ebj_without_direct_2.jld2"; grid_dimension = 3, grid_refinement_max = 2, solver_config = :gmres_bicgstab_ebj_without_direct)

    ############################################################################
    # Runtime summary
    ############################################################################

    println("\nMOSFET benchmark runtimes:")
    println("tdirect_all_2d                            = $tdirect_all_2d seconds")
    println("tgmres_ilu_part_2d                        = $tgmres_ilu_part_2d seconds")
    println("tgmres_amg_part_2d                        = $tgmres_amg_part_2d seconds")
    println("tbicgstab_ilu_part_2d                     = $tbicgstab_ilu_part_2d seconds")
    println("tbicgstab_amg_part_2d                     = $tbicgstab_amg_part_2d seconds")
    println("tgmres_ebj_2d                             = $tgmres_ebj_2d seconds")
    println("tbicgstab_ebj_2d                          = $tbicgstab_ebj_2d seconds")
    println("tgmres_bicgstab_ebj_2d                    = $tgmres_bicgstab_ebj_2d seconds")
    println("tdirect_pardiso_3d                        = $tdirect_pardiso_3d seconds")
    println("tgmres_bicgstab_ebj_3d                    = $tgmres_bicgstab_ebj_3d seconds")
    println("tgmres_bicgstab_ebj_without_direct_3d     = $tgmres_bicgstab_ebj_without_direct_3d seconds")

    return nothing
end # several_benchmarks

end # module
