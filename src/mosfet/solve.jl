################################################################################
# Solve system in equilibrium
################################################################################
# Taken from ChargeTransport.jl and adapted (charge_neutrality and lambda_embedding as initial value options added, benchmark_setting added)
function equilibrium_solve!(ctsys::System, initial_value, grid_dimension; inival = VoronoiFVM.unknowns(ctsys.fvmsys, inival = 0.0), control = VoronoiFVM.NewtonControl(), benchmark_data = BenchmarkSeries(), vacancyEnergyCalculation::Bool = false, verbose::Bool = false, yabstol::Float64 = 1.0e-2, ytol::Float64 = 1.0e-4, maxiter::Int64 = 15) # last three are extended-only keywords for vacancyEnergyCalculation = true

    # by default vacancyEnergyCalculation false
    return _equilibrium_solve!(Val(vacancyEnergyCalculation), ctsys, initial_value, grid_dimension; inival = inival, control = control, benchmark_data = benchmark_data, verbose = verbose, yabstol = yabstol, ytol = ytol, maxiter = maxiter)

end

# Base implementation of equilibrium_solve: vacancyEnergyCalculation = false
function _equilibrium_solve!(::Val{false}, ctsys::System, initial_value, grid_dimension; inival, control, benchmark_data, verbose = verbose, yabstol, ytol, maxiter)

    ctsys.fvmsys.physics.data.calculationType = InEquilibrium
    grid = ctsys.fvmsys.grid

    data = ctsys.fvmsys.physics.data
    params = ctsys.fvmsys.physics.data.params
    paramsnodal = ctsys.fvmsys.physics.data.paramsnodal
    bnode = grid[BFaceNodes]
    ipsi = data.index_psi
    iphin = 1
    iphip = 2
    (; k_B, q) = data.constants

    # We set zero voltage for each charge carrier at all outer boundaries for equilibrium calculations.
    for ibreg in grid[BFaceRegions]
        set_contact!(ctsys, ibreg, Δu = 0.0)
    end

    sol = inival

    # alternative initial value calculation: use charge neutrality condition
    if initial_value == :charge_neutrality

        U_T = k_B * params.temperature / q

        # note: everywhere the same, one material
        Nc = params.densityOfStates[iphin, 1]
        Nv = params.densityOfStates[iphip, 1]
        Ec = params.bandEdgeEnergy[iphin, 1]
        Ev = params.bandEdgeEnergy[iphip, 1]

        n_intr = sqrt(Nc * Nv * exp(-(Ec - Ev) / (k_B * params.temperature)))

        nnodes = num_nodes(grid)
        coord = grid[Coordinates]
        @local_unitfactors μm

        for i in 1:nnodes
            # doping values defined cell-wise rather than at the grid nodes, where the electrostatic potential is defined
            # doping region is determined from the grid coordinates, separately for the 2D and 3D grids
            x = coord[1, i]
            y = coord[2, i]

            # region boundaries hard-coded for the specific grids used here, could be generalized in the future!!
            if grid_dimension == 2
                region = (((-1.5μm <= x <= -0.8μm) || (0.8μm <= x <= 1.5μm)) && (-0.2μm <= y <= 0.0μm)) ? 1 : 2

            elseif grid_dimension == 3
                z = coord[3, i]
                region = (((-1.5μm <= x <= -0.8μm) || (0.8μm <= x <= 1.5μm)) && (0μm <= y <= 7.5μm) && (-0.2μm <= z <= 0.0μm)) ? 1 : 2

            end

            # C is a scalar at each node, but evaluating in loop yields the spatial doping profile
            C = params.doping[iphin, region] - params.doping[iphip, region]

            inival[ipsi, i] = (Ec + Ev) / (2q) - 0.5 * U_T * log(Nc / Nv) + U_T * asinh(C / (2 * n_intr))
        end

        sol = ChargeTransport.solve(ctsys; inival = inival, control = control)

        history = VoronoiFVM.history(sol)

        push!(
            benchmark_data.equilibrium,
            BenchmarkStep(
                value = 1.0,
                newton_iters = length(history.updatenorm),
                tlinsolve = history.tlinsolve,
                tlinsolve_setup = history.tlinsolve_setup,
                tlinsolve_solve = history.tlinsolve_solve
            )
        )

    elseif initial_value == :lambda_embedding
        # lambda1: 0 (linear problem) to 1 (full nonlinear problem)
        range_lambda1 = [0.0, 1.0]

        # callback executed before each embedding step
        control.pre = (solution, lambda1) -> begin
            if control.verbose == "n"
                println("λ1 = $(lambda1)")
            end
            ctsys.fvmsys.physics.data.λ1 = lambda1
            nothing
        end

        control.post = (solution, oldsolution, lambda1, Δp) -> begin
            history = VoronoiFVM.history(solution)
            push!(
                benchmark_data.equilibrium, BenchmarkStep(
                    value = lambda1,
                    newton_iters = length(history.updatenorm),
                    tlinsolve = history.tlinsolve,
                    tlinsolve_setup = history.tlinsolve_setup,
                    tlinsolve_solve = history.tlinsolve_solve
                )
            )
            nothing
        end

        # solve the stationary problem with parameter embedding (bias continuation)
        sol_allLambdaValues = ChargeTransport.solve(ctsys; inival = inival, control = control, embed = range_lambda1)
        sol = sol_allLambdaValues.u[end]

    else
        error("Unknown initial_value = $initial_value")
    end

    for ibreg in grid[BFaceRegions]
        # here we assume that in multidimensions, we receive a constant value of the electric potential at the boundary
        # check for applications, where this is not the case
        bψVal = view(sol[ipsi, :], subgrid(grid, [ibreg], boundary = true))[1]
        params.bψEQ[ibreg] = bψVal
    end

    # calculate equilibrium densities (especially needed for Schottky boundary conditions)
    for icc in data.electricCarrierList
        for ibreg in grid[BFaceRegions]
            Ncc = params.bDensityOfStates[icc, ibreg] + paramsnodal.densityOfStates[icc, bnode[ibreg]]
            Ecc = params.bBandEdgeEnergy[icc, ibreg] + paramsnodal.bandEdgeEnergy[icc, bnode[ibreg]]

            eta = params.chargeNumbers[icc] / (k_B * params.temperature / q) * ((sol[icc, bnode[ibreg]] - sol[ipsi, bnode[ibreg]]) + Ecc / q)
            params.bDensityEQ[icc, ibreg] = Ncc * data.F[icc](eta)
        end
    end

    # set now calculationType to outOfEquilibrium for further calculations
    data.calculationType = OutOfEquilibrium

    ## put here back the homogeneous Neumann boundary conditions.
    for iicc in data.ionicCarrierList

        icc = iicc.ionicCarrier # species number chosen by user

        for ibreg in 1:grid[NumBFaceRegions]
            ctsys.fvmsys.boundary_factors[icc, ibreg] = 0.0
            ctsys.fvmsys.boundary_values[icc, ibreg] = 0.0
        end

    end
    for iicc in data.trapCarrierList

        icc = iicc.trapCarrier # species number chosen by user

        for ibreg in 1:grid[NumBFaceRegions]
            ctsys.fvmsys.boundary_factors[icc, ibreg] = 0.0
            ctsys.fvmsys.boundary_values[icc, ibreg] = 0.0
        end

    end

    ## this flag is needed, when the user e.g. defines customized contactVoltageFunctions
    data.generationComplete = true

    # save changes on fvmsys of VoronoiFVM likewise in ctsys.data
    ctsys.data = ctsys.fvmsys.physics.data

    return sol

end


################################################################################
# Solve system for applied gate voltage (V_G: 0V to 5V, V_D = 0V)
################################################################################
function gate_solve!(ctsys::System, biasValue_gate_range; inival = VoronoiFVM.unknowns(ctsys.fvmsys, inival = 0.0), control = VoronoiFVM.NewtonControl(), p = ChargeTransport.Params(), benchmark_data = BenchmarkSeries())

    # callback executed before each embedding step
    control.pre = (solution, biasValue_gate) -> begin
        println("bias value at gate: Δu = ", biasValue_gate, "V")
        set_contact!(ctsys, p.bregion_gate, Δu = biasValue_gate)
        nothing
    end

    control.post = (solution, oldsolution, biasValue_gate, Δp) -> begin
        history = VoronoiFVM.history(solution)
        push!(
            benchmark_data.gate, BenchmarkStep(
                value = biasValue_gate,
                newton_iters = length(history.updatenorm),
                tlinsolve = history.tlinsolve,
                tlinsolve_setup = history.tlinsolve_setup,
                tlinsolve_solve = history.tlinsolve_solve
            )
        )
        nothing
    end

    # solve the stationary problem with parameter embedding (bias continuation)
    solution_allBiasValues = ChargeTransport.solve(ctsys; inival = inival, control = control, embed = biasValue_gate_range)
    solution = solution_allBiasValues.u[end]

    return solution
end

################################################################################
# Solve system for applied drain voltage (V_G = 5V, V_D: 0V to 5V)
################################################################################
function drain_solve!(ctsys::System, biasValue_drain_range, grid_dimension; inival = VoronoiFVM.unknowns(ctsys.fvmsys, inival = 0.0), control = VoronoiFVM.NewtonControl(), p = ChargeTransport.Params(), benchmark_data = BenchmarkSeries())

    # prepare getting IV data
    IV_drain = Float64[]
    biasValues_drain = Float64[]
    factory = VoronoiFVM.TestFunctionFactory(ctsys.fvmsys)
    tf = VoronoiFVM.testfunction(factory, [p.bregion_drain], [p.bregion_source])

    # callback executed before each embedding step
    control.pre = (solution, biasValue_drain) -> begin
        println("bias value at drain: Δu = ", biasValue_drain, "V")
        set_contact!(ctsys, p.bregion_drain, Δu = biasValue_drain)
        nothing
    end

    # callback executed after each embedding step
    control.post = (solution, oldsolution, biasValue_drain, Δp) -> begin
        # get IV data
        IEdge = VoronoiFVM.integrate_∇TxFlux(ctsys.fvmsys, tf, solution)
        current = IEdge[p.iphin] + IEdge[p.iphip]
        if grid_dimension == 2
            push!(IV_drain, abs(p.zaus * current)) # in 2D, we need to multiply with the depth of the device (here: 1 m) to get a current value
        elseif grid_dimension == 3
            push!(IV_drain, abs(2 * current)) # in 3D, we need to multiply with 2, because we only simulate half of the device (see generate_3Dgrid_half)
        end
        push!(biasValues_drain, biasValue_drain)

        # get runtimes linear Solver
        history = VoronoiFVM.history(solution)
        push!(
            benchmark_data.drain, BenchmarkStep(
                value = biasValue_drain,
                newton_iters = length(history.updatenorm),
                tlinsolve = history.tlinsolve,
                tlinsolve_setup = history.tlinsolve_setup,
                tlinsolve_solve = history.tlinsolve_solve
            )
        )
        nothing
    end

    # solve the stationary problem with parameter embedding (bias continuation)
    solution_allBiasValues = ChargeTransport.solve(ctsys; inival = inival, control = control, embed = biasValue_drain_range)
    solution = solution_allBiasValues.u[end]

    return solution, biasValues_drain, IV_drain
end
