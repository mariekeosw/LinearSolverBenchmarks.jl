### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# ╔═╡ c31f4cc9-84ed-4897-8678-133ca8fc3881
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
end

# ╔═╡ f38de270-0731-11f1-b853-e7a961d52743
begin
    # Grid, Parameters, ...
    using LinearSolverBenchmarks: Mosfet
    using ExtendableGrids
    using ChargeTransport
    using VoronoiFVM # for current density
    using LessUnitful: @unitfactors

    # Visualization
    using GridVisualize
    using PlutoVista
    using CairoMakie
    using LaTeXStrings

    # Data Import
    using JLD2

    # Table of Contents
    using PlutoUI
end

# ╔═╡ e2cbe387-93bb-4974-aef0-069e2887b654
md"""
# 3D Result Plots
"""

# ╔═╡ 96ee39ad-8e47-48ff-a72d-d33c89f8e1a5
TableOfContents()

# ╔═╡ cc8a2d3d-4845-4703-91bd-c92f0ee95414
md"""
## Grid Plot
- shorter gate electrode for 3D effect
- half grid, to shorten the simulation time (current is in the end multiplied by 2)
- refinement as in 2D
"""

# ╔═╡ 7c4f4d1c-5138-45c8-9690-1be3410105aa
begin
    @load "data_mosfet_3d.jld2" grid_refinement solution_eq solution biasValues_drain IV_drain

    grid = Mosfet.generate_3Dgrid_half(grid_refinement)

    # For better readability
    μm = 1.0e-6
    grid_um = deepcopy(grid)
    grid_um[Coordinates] ./= μm

    GridVisualize.gridplot(
        grid_um,
        Plotter = CairoMakie,
        xlabel = "width [1e-6 m]",
        ylabel = "length [1e-6 m]",
        zlabel = "height [1e-6 m]",
		# Ansicht Gitter
		resolution = (800, 500),
        yplanes = [7.5],
        azim = -60, # Drehung, Default -60
        elev = 30, # Default
        outlinealpha = 0.4,
		fontsize = 50
    )
end

# ╔═╡ 5e50ae16-7b3e-4758-be81-d01e9577541e
let
	# Doping - but not working nicely with CairoMakie
	@load "data_mosfet_3d.jld2" grid_refinement solution_eq solution biasValues_drain IV_drain

    grid = Mosfet.generate_3Dgrid_half(grid_refinement)

    # For better readability
    μm = 1.0e-6
    grid_um = deepcopy(grid)
    grid_um[Coordinates] ./= μm

    GridVisualize.gridplot(
        grid_um,
        Plotter = CairoMakie,
        xlabel = "width [µm]",
        ylabel = "length [µm]",
        zlabel = "height [µm]",
		resolution = (1000, 500),
        yplanes = [5.0],
        azim = 90, # Drehung, Default -60
        elev = 30, # Default
        outlinealpha = 0.0,
		linewidth = 0,
		planealpha = 1.0,
	)
end

# ╔═╡ e6602e2f-4050-4d88-aeb3-38db2636b751
md"""
## Solution Plots in Equilibrium

### Electrostatic Potential

"""

# ╔═╡ 2913e714-450f-48ac-81c2-ec4b41068ea4
let
    # Electrostatic Potential Equilibrium
    @load "data_mosfet_3d.jld2" grid_refinement solution_eq solution biasValues_drain IV_drain

    psi = Vector(solution_eq[3, :])

    scalarplot(
        grid_um,
        psi;
        Plotter = CairoMakie,
        resolution = (1000, 500),
        xlabel = "width [µm]",
        ylabel = "length [µm]",
        zlabel = "height [µm]",
        colorbarlabel = "potential [V]",
        levels = 10,
        yplanes = [7.49], # Es muss ein GItterpunkt getroffen werden
        azim = -60, # 90, 0 - dann hat man eine frontale Sicht
        elev = 30, # Default
        colorbar = :vertical,
		fontsize = 40
    )
end

# ╔═╡ 2913c5d1-6043-40b2-9096-bb6cf6a39e5c
md"""
### Electron density
"""

# ╔═╡ 88e5ed8c-824d-45b7-ae86-62779f77fbf3
let
    # Electron Density Equilibrium
    @load "data_mosfet_3d.jld2" grid_refinement solution_eq solution biasValues_drain IV_drain

    # Get MOSFET parameter setup and constants
    parameter_set = Mosfet.tesca_setup
    p = parameter_set()
    constants = ChargeTransport.teSCA_constants

    # Calculate electron density
    nn = p.Nc .* exp.(p.zn * (constants.q * (solution_eq[p.iphin, :] .- solution_eq[p.ipsi, :]) .+ p.Ec) ./ (constants.k_B * p.T))
    nn_cm = nn * 1.0e-6

    # log scaling of colormap as in 2D (Electron Density)
    ticks = [1.0e7, 1.0e9, 1.0e11, 1.0e13, 1.0e15, 1.0e17, 1.0e19]

    scalarplot(
        grid_um,
        nn_cm;
        Plotter = CairoMakie,
        resolution = (1000, 500),
        xlabel = "width [µm]",
        ylabel = "length [µm]",
        zlabel = "height [µm]",
        colorscale = :log,
        colorbarticks = ticks,
        colorbarlabel = "density [cm⁻³]",
        yplanes = [7.49], # Es muss ein GItterpunkt getroffen werden
        azim = -60, # Drehung, Default -60
        elev = 30, # Default
        colorbar = :vertical,
		fontsize = 40
    )

end

# ╔═╡ 1ac0ec0b-f76f-463c-a68a-17b839b5fffc
md"""
### Hole Denisty
"""

# ╔═╡ 7b1fc13b-df18-44d0-9da3-76c10d1ec788
let
    # Hole Density Equilibrium
    @load "data_mosfet_3d.jld2" grid_refinement solution_eq solution biasValues_drain IV_drain

    # Get MOSFET parameter setup and constants
    parameter_set = Mosfet.tesca_setup
    p = parameter_set()
    constants = ChargeTransport.teSCA_constants

    # Calculate hole density
    np = p.Nv .* exp.(p.zp * (constants.q * (solution_eq[p.iphip, :] .- solution_eq[p.ipsi, :]) .+ p.Ev) ./ (constants.k_B * p.T))
    np_cm = np * 1.0e-6

    # log scaling of colormap as in 2D (Hole Density)
    ticks = [1.0e4, 1.0e6, 1.0e8, 1.0e10, 1.0e12, 1.0e14, 1.0e16]

    scalarplot(
        grid_um,
        np_cm;
        Plotter = CairoMakie,
        resolution = (1000, 500),
        xlabel = "width (µm)",
        ylabel = "length (µm)",
        zlabel = "height (µm)",
        colorscale = :log,
        colorbarticks = ticks,
        colorbarlabel = "density [cm⁻³]",
        levels = 10,
        yplanes = [7.49], # Es muss ein GItterpunkt getroffen werden
        azim = -60, # Drehung, Default -60
        elev = 30, # Default
        colorbar = :vertical,
		fontsize = 40
    )
end

# ╔═╡ 65b1550f-918b-46bf-a133-9e988a721f07
md"""
## Solution Plots for Full Solution
### Electrostatic Potential
"""

# ╔═╡ e680170a-4d40-47b6-94ba-36d03849a70e
let
    # Electrostatic Potential
    @load "data_mosfet_3d.jld2" grid_refinement solution_eq solution biasValues_drain IV_drain

    psi = Vector(solution[3, :])

    scalarplot(
        grid_um,
        psi;
        Plotter = CairoMakie,
        resolution = (600, 500),
        xlabel = "width [1e-6 m]",
        ylabel = "length [1e-6 m]",
        zlabel = "height [1e-6 m]",
        colorbarlabel = "potential [V]",
        levels = 10,
        yplanes = [7.49], # Es muss ein GItterpunkt getroffen werden
        azim = -60, # 90, 0 - dann hat man eine frontale Sicht
        elev = 30, # Default
        colorbar = :vertical,
		fontsize = 50
    )
end

# ╔═╡ 4dab8657-d75e-4c19-8d94-bd5318520273
md"""
### Electron density
"""

# ╔═╡ 4a043fe4-d5c2-43ad-8c01-eebd3e29e36d
let
    # Electron Density
    @load "data_mosfet_3d.jld2" grid_refinement solution_eq solution biasValues_drain IV_drain

    # Get MOSFET parameter setup and constants
    parameter_set = Mosfet.tesca_setup
    p = parameter_set()
    constants = ChargeTransport.teSCA_constants

    # Calculate electron density
    nn = p.Nc .* exp.(p.zn * (constants.q * (solution[p.iphin, :] .- solution[p.ipsi, :]) .+ p.Ec) ./ (constants.k_B * p.T))
    nn_cm = nn * 1.0e-6

    # log scaling of colormap as in 2D (Electron Density)
    ticks = [1.0e7, 1.0e9, 1.0e11, 1.0e13, 1.0e15, 1.0e17, 1.0e19]

    scalarplot(
        grid_um,
        nn_cm;
        Plotter = CairoMakie,
        resolution = (600, 500),
        xlabel = "width [1e-6 m]",
        ylabel = "length [1e-6 m]",
        zlabel = "height [1e-6 m]",
        colorscale = :log,
        colorbarticks = ticks,
        colorbarlabel = "density [cm⁻³]",
        yplanes = [7.49], # Es muss ein GItterpunkt getroffen werden
        azim = -60, # Drehung, Default -60
        elev = 30, # Default
        colorbar = :vertical,
		fontsize = 50
    )
end

# ╔═╡ 1e62f4b2-1c7e-4177-acaf-66ec5a323567
md"""
### Hole density
"""

# ╔═╡ 81092801-4f4f-4a73-8348-1978c0440482
let
    # Hole Density Equilibrium
    @load "data_mosfet_3d.jld2" grid_refinement solution_eq solution biasValues_drain IV_drain

    # Get MOSFET parameter setup and constants
    parameter_set = Mosfet.tesca_setup
    p = parameter_set()
    constants = ChargeTransport.teSCA_constants

    # Calculate hole density
    np = p.Nv .* exp.(p.zp * (constants.q * (solution[p.iphip, :] .- solution[p.ipsi, :]) .+ p.Ev) ./ (constants.k_B * p.T))
    np_cm = np * 1.0e-6

    # log scaling of colormap as in 2D (Hole Density)
    ticks = [1.0e4, 1.0e6, 1.0e8, 1.0e10, 1.0e12, 1.0e14, 1.0e16]

    scalarplot(
        grid_um,
        np_cm;
        Plotter = CairoMakie,
        resolution = (800, 500),
        xlabel = "width [µm]",
        ylabel = "length [µm]",
        zlabel = "height [µm]",
        colorscale = :log,
        colorbar = :vertical,
        colorbarticks = ticks,
        colorbarlabel = "density [cm⁻³]",
        levels = 10,
        yplanes = [7.49], # Es muss ein GItterpunkt getroffen werden
        azim = -60, # Drehung, Default -60
        elev = 30, # Default
		fontsize = 35

    )
end

# ╔═╡ fc9ff41e-3d32-4568-a04c-aedfa8d45d46
md"""
### Current Density
"""

# ╔═╡ a91eb036-cc3a-4c48-aef1-cb2f9e92bedf
let
    # Current Density
    @load "data_mosfet_3d.jld2" grid_refinement solution_eq solution biasValues_drain IV_drain

    #############################################################################
    # Define VoronoiFVM System
    #############################################################################
    parameter_set = tesca_setup
    p = parameter_set()
    constants = ChargeTransport.teSCA_constants

    # Initialize Data instance
    data = Data(grid, p.numberOfCarriers, constants = constants)

    # model
    data.modelType = Stationary

    # statistics
    data.F .= Boltzmann # Tesca Manual S.84 IFERMI = 0

    # recombination
    data.bulkRecombination = set_bulk_recombination(;
        iphin = p.iphin, iphip = p.iphip,
        bulk_recomb_Auger = true,
        bulk_recomb_radiative = false,
        bulk_recomb_SRH = true,
    )

    # boundary model
    data.boundaryType[p.bregion_gate] = GateContact
    data.boundaryType[p.bregion_drain] = OhmicContact
    data.boundaryType[p.bregion_source] = OhmicContact
    data.boundaryType[p.bregion_bulk] = OhmicContact

    # flux discretization - depends on statistic (Boltzmann -> Scharfetter Gummel)
    data.fluxApproximation .= ScharfetterGummel

    # add parameters to data
    data.params = ChargeTransport.Params(p)

    # Definition ChargeTransport System
    ctsys = System(grid, data, unknown_storage = :sparse)
    #############################################################################

    nf = -VoronoiFVM.nodeflux(ctsys.fvmsys, solution)

    j = nf[:, p.iphin, :]

    # Eingabewerte von 2D
    vmin = 1.0e-8
    vmax = 1.0e4

    jx = Vector(j[1, :])
    jy = Vector(j[2, :])
    jz = Vector(j[3, :])

    jabs = sqrt.(jx .^ 2 .+ jy .^ 2 .+ jz .^ 2)

    levels = 10.0 .^ range(log10(vmin), log10(vmax), length = 25)

    scalarplot(
        grid_um,
        jabs;
        Plotter = CairoMakie,
        resolution = (700, 500),
        xlabel = "width [µm]",
        ylabel = "length [µm]",
        zlabel = "height [µm]",
        colorscale = :log,
        colorbar = :vertical,
        colorbarticks = 10.0 .^ (-8:2:4),
        colorbarlabel = "|j| [A/cm²]",
        levels = levels,
        yplanes = [7.49],
        azim = -60,
        elev = 30
    )

end

# ╔═╡ 193e33d5-ce41-4b91-8e48-2ea1230c1122
md"""
### IV Curve
"""

# ╔═╡ e783f0e7-cbec-4117-b0bb-6aa2f89257d9
let
    # IV Curve
    @load "data_mosfet_3d.jld2" grid_refinement solution_eq solution biasValues_drain IV_drain

    scalarplot(
        biasValues_drain,
        IV_drain;
        Plotter = CairoMakie,
        resolution = (600, 500),
        xlabel = "voltage drain [V]",
        ylabel = "total current drain [A]",
		xticks = 0:1:5,
		title = " ",
		fontsize = 60,
		markershape = :circle,
		markersize = 14,
		markevery = 1, # every measurement is shown
		color = "#1f77b4", # matplotlib blue, analog to 2d plots
    )
end

# ╔═╡ Cell order:
# ╟─e2cbe387-93bb-4974-aef0-069e2887b654
# ╟─96ee39ad-8e47-48ff-a72d-d33c89f8e1a5
# ╠═c31f4cc9-84ed-4897-8678-133ca8fc3881
# ╠═f38de270-0731-11f1-b853-e7a961d52743
# ╟─cc8a2d3d-4845-4703-91bd-c92f0ee95414
# ╠═7c4f4d1c-5138-45c8-9690-1be3410105aa
# ╠═5e50ae16-7b3e-4758-be81-d01e9577541e
# ╟─e6602e2f-4050-4d88-aeb3-38db2636b751
# ╠═2913e714-450f-48ac-81c2-ec4b41068ea4
# ╟─2913c5d1-6043-40b2-9096-bb6cf6a39e5c
# ╠═88e5ed8c-824d-45b7-ae86-62779f77fbf3
# ╟─1ac0ec0b-f76f-463c-a68a-17b839b5fffc
# ╠═7b1fc13b-df18-44d0-9da3-76c10d1ec788
# ╟─65b1550f-918b-46bf-a133-9e988a721f07
# ╠═e680170a-4d40-47b6-94ba-36d03849a70e
# ╟─4dab8657-d75e-4c19-8d94-bd5318520273
# ╠═4a043fe4-d5c2-43ad-8c01-eebd3e29e36d
# ╟─1e62f4b2-1c7e-4177-acaf-66ec5a323567
# ╠═81092801-4f4f-4a73-8348-1978c0440482
# ╟─fc9ff41e-3d32-4568-a04c-aedfa8d45d46
# ╟─a91eb036-cc3a-4c48-aef1-cb2f9e92bedf
# ╟─193e33d5-ce41-4b91-8e48-2ea1230c1122
# ╠═e783f0e7-cbec-4117-b0bb-6aa2f89257d9
