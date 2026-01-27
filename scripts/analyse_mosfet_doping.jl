#=
This file complements "script/run_mosfet_benchmark.jl". While the main
benchmark script uses region-wise doping, this script replaces the region-wise
doping values by a node-wise doping profile.

The node-wise doping profile is imported from TeSCA data and is therefore only
tested for the two-dimensional setup. In this setting, the grid has to match the
corresponding TeSCA grid (default ref=0).

For this node-wise doping configuration, the charge-neutrality initial value
does not converge reliably. Therefore, lambda embedding is used as the preferred
initialization strategy.
=#

module MosfetDopingAnalysis

using ChargeTransport
using VoronoiFVM
using ExtendableGrids
using LinearSolverBenchmarks: Mosfet
using JLD2

# How to use this function?
# 1) includet("scripts/analyse_mosfet_doping.jl")
# 2) MosfetDopingAnalysis.main_nodewise("post_process/figures/mosfet_2d_tesca/tesca_files/data_mosfet_ref0_doping_dio.dat", 0)
function main_nodewise(
        tesca_doping_file::String,
        grid_refinement;
        ###############################
        Documentation = true,
        ###############################
        grid_dimension = 2,                # fixed
        ###############################
        initial_value = :lambda_embedding, # fixed
        strategy = nothing,
        ###############################
        parameter_set = Mosfet.tesca_setup
    )

    # Initialize storage for benchmark data (runtimes, iterations)
    benchmark_data = Mosfet.BenchmarkSeries()

    ################################################################################
    if Documentation
        println("Set up grid and regions")
    end
    ################################################################################

    # Select the two-dimensional grid.
    # The grid must be compatible with the imported TeSCA doping profile.
    grid = Mosfet.generate_2Dgrid(grid_refinement)

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

    # Import parameters from parameter file (TeSCA setup)
    p = parameter_set()

    # Initialize Data instance
    data = Data(grid, p.numberOfCarriers, constants = constants)

    data.modelType = Stationary

    data.F .= Boltzmann # Tesca Manual S.84 IFERMI = 0

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

    # Add standard parameter set to the simulation data
    data.params = ChargeTransport.Params(p)

    # Disable region-wise doping by resetting all doping values to zero
    data.params.doping[p.iphin, p.region_n] = 0
    data.params.doping[p.iphip, p.region_p] = 0

    # Load node-wise doping profile from a TeSCA reference simulation with EPSSi ≈ 0.
    # In this limit, Poisson's equation reduces to
    #     0 = n - p + C
    # so the doping profile can be extracted as
    #     C = -(n - p).
    data.paramsnodal = ChargeTransport.ParamsNodal(grid, p, tesca_doping_file)

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

    # Save simulation results for post-processing and visualization
    @save "data_mosfet_nodal_doping.jld2" grid_refinement solution_eq solution biasValues_drain IV_drain

    return num_nodes, testval, testval_phy, benchmark_summary.newton_iters, benchmark_summary.tlinsolve, benchmark_summary.tlinsolve_setup, benchmark_summary.tlinsolve_solve
end # function main

end # module
