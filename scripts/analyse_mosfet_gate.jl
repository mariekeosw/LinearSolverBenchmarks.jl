#=
This file complements "script/run_mosfet_benchmark.jl" and generates data for
the gate-voltage-dependent output characteristics of the MOSFET.

The simulation setup is the same as in the main benchmark script. However, the
simulation is repeated for different maximum gate voltages.

In addition, the drain-voltage continuation is adapted to store the current
values at every embedding step. These values are used to construct the
current-voltage curves for the different gate voltages.
=#
module MosfetGateAnalysis

using ChargeTransport
using VoronoiFVM          # for IV Curve, Current Density Plot
using ExtendableGrids
using LinearSolverBenchmarks: Mosfet
using JLD2

# How to use this function?
# 1) includet("scripts/analyse_mosfet_gate.jl")
# 2) MosfetGateAnalysis.main()
function main(;
        Documentation = true,
        ###############################
        grid_dimension = 2,
        grid_refinement = 0,  # 0 with xx nodes, 1 with xx nodes, ...
        ###############################
        initial_value = :charge_neutrality, # initial value for the equilibrium run (either: :charge_neutrality or :lambda_embedding)
        strategy = nothing, # Default for sparse matrices is UMFPACKFactorization()
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
    control.log = true # Record history, ist per Default auf false
    # Newton control parameters
    control.maxiters = 20
    control.abstol = 1.0e-10 # 1.0e-10
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

    # Different max biases, each bias will procuce one line in the final plot
    maxbias_gate = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]

    biasValues_drain = Dict{Float64, Vector{Float64}}()
    IV_drain = Dict{Float64, Vector{Float64}}()


    for gate_bias in maxbias_gate

        # Start each gate curve from equilibrium
        inival_gate = copy(solution_eq)

        # Apply fixed gate voltage
        biasValue_gate_range = [0.0, gate_bias]
        solution_gate = Mosfet.gate_solve!(ctsys, biasValue_gate_range, inival = inival_gate, control = control, p = p, benchmark_data = benchmark_data)

        # drain bias: start at 0 V and increase up to 5 V
        inival_drain = copy(solution_gate)
        biasValue_drain_range = [0.0, 5.0]
        solution_drain, drain_biases, drain_currents = Mosfet.drain_solve!(ctsys, biasValue_drain_range, grid_dimension, inival = inival_drain, control = control, p = p, benchmark_data = benchmark_data)

        biasValues_drain[gate_bias] = drain_biases
        IV_drain[gate_bias] = drain_currents

        # Reset contacts to equilibrium biases
        set_contact!(ctsys, p.bregion_gate, Δu = 0.0)
        set_contact!(ctsys, p.bregion_drain, Δu = 0.0)
    end
    if Documentation
        println("*** done\n")
    end

    @save "data_mosfet_2d_gate_property.jld2" maxbias_gate biasValues_drain IV_drain

    return nothing
end # function main

end # module
