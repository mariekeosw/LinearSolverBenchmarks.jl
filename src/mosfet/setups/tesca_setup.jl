# Parameters are taken from the Tesca Manual

@kwdef struct tesca_setup

    # constants
    constants = ChargeTransport.teSCA_constants

    # used unit factors
    K = ufac"K"
    m = ufac"m"
    cm = ufac"cm"
    V = ufac"V"
    s = ufac"s"

    eV = constants.q * V

    ########## charge carriers ##########
    iphin = 1                            # quasi Fermi potential for electrons
    iphip = 2                            # quasi Fermi potential for holes
    ipsi = 3                             # electrostatic potential
    numberOfCarriers = 2

    ########## device geometry ##########
    # region numbers
    region_n = 1 # n-doped region
    region_p = 2 # p-doped region
    numberOfRegions = 2

    # boundary region numbers
    bregion_gate = 1
    bregion_drain = 2
    bregion_source = 3
    bregion_bulk = 4
    bregion_neumann = 5
    numberOfBoundaryRegions = 5

    zaus = 15.0e-4 * cm                                 # depth of the device [cm], see mosS.dio file
    thickness_ox = 0.044e-4 * cm                        # oxide thickness on gate [cm], see mosS.dio file

    ########## physical values of Si at room temperature ##########

    # charge numbers
    zn = -1
    zp = 1

    # temperature
    T = 300.0 * K

    # band edge energies
    Ec = 0.562 * eV                                     # conduction band-edge energy, Tesca Manual S. 83
    Ev = -0.562 * eV                                    # valence band-edge energy, Tesca Manual S. 83

    # density of states
    Nc = 2.86e19 / (cm^3)                               # conduction band density of states, Tesca Manual S.85
    Nv = 3.1e19 / (cm^3)                                # valence band densitiy of states, Tesca Manual S.85
    Nss = 3.2251421571687153e11 / (cm^2)                # surface state density

    # density of states
    mun = 1030.0 * (cm^2) / (V * s)                     # electron mobility, ioffe \leq 1400, Tesca Manual S.92: 1030
    mup = 495.0 * (cm^2) / (V * s)                      # hole mobility \ioffe \leq 450, Tesca Manual S. 92: 495

    # relative dielectric permittivities
    εr = 11.67 * 1.0                                    # relative dielectric permittivity of Si, Tesca Manual S.21: 11.67
    εr_ox = 3.8 * 1.0                                   # relative dielectric permittivity of SiO2 at gate, Tesca Manual S. 21. 3.8

    # Recombination parameters (from Tesca)
    # Radiative = 0 # Tesca Manual S.127
    Auger_n = 2.8e-31                         # Tesca Manual S. 127: 2.8d-31
    Auger_p = 9.9e-32                         # Tesca Manual S. 127: 9.9d-32
    SRH_TrapDensity_n = 1.09e10 / cm^3        # Tesca Manual S. 128: 1.09d10 1/cm³
    SRH_TrapDensity_p = 1.09e10 / cm^3        # Tesca Manual S. 128: 1.09d10 1/cm³
    SRH_LifeTime_n = 2.0e-4 * s               # Tesca Manual S. 128, tau_n0 = 2e-4 s
    SRH_LifeTime_p = 2.0e-6 * s               # Tesca Manual S. 128, trau_p0 = 2e-6 s
    SRH_Velocity_n = 5.0 * (cm / s)           # Tesca Manual S. 129, VREN = 5.d0 cm/s
    SRH_Velocity_p = 5.0 * (cm / s)           # Tesca Manual S. 129, VREP = 5.d0 cm/s

    # Doping
    Na = 1.0e16 / cm^3
    Nd = 1.0e20 / cm^3
end

# Now we create an ChargeTransport.Params object
function ChargeTransport.Params(p::tesca_setup)

    params = ChargeTransport.Params(
        p.numberOfRegions,
        p.numberOfBoundaryRegions,
        p.numberOfCarriers
    )

    params.temperature = p.T
    params.chargeNumbers[p.iphin] = p.zn
    params.chargeNumbers[p.iphip] = p.zp

    for ireg in 1:p.numberOfRegions # region data
        params.dielectricConstant[ireg] = p.εr * p.constants.ε_0

        # effective DOS, band-edge energy and mobilities
        params.densityOfStates[p.iphin, ireg] = p.Nc
        params.densityOfStates[p.iphip, ireg] = p.Nv
        params.bandEdgeEnergy[p.iphin, ireg] = p.Ec
        params.bandEdgeEnergy[p.iphip, ireg] = p.Ev
        params.mobility[p.iphin, ireg] = p.mun
        params.mobility[p.iphip, ireg] = p.mup

        params.recombinationSRHLifetime[p.iphin, ireg] = p.SRH_LifeTime_n
        params.recombinationSRHLifetime[p.iphip, ireg] = p.SRH_LifeTime_p
        params.recombinationSRHTrapDensity[p.iphin, ireg] = p.SRH_TrapDensity_n
        params.recombinationSRHTrapDensity[p.iphip, ireg] = p.SRH_TrapDensity_p
        params.recombinationAuger[p.iphin, ireg] = p.Auger_n
        params.recombinationAuger[p.iphip, ireg] = p.Auger_p
    end

    #gate-specific recombination
    params.recombinationSRHvelocity[p.iphin, p.bregion_gate] = p.SRH_Velocity_n
    params.recombinationSRHvelocity[p.iphip, p.bregion_gate] = p.SRH_Velocity_p

    # gate-specific parameters
    params.dielectricConstantOxideGate[p.bregion_gate] = p.εr_ox * p.constants.ε_0
    params.thicknessOxideGate[p.bregion_gate] = p.thickness_ox
    params.surfaceChargeDensityGate[p.bregion_gate] = p.constants.q * p.Nss

    # doping                                                   #     n    p
    params.doping[p.iphin, p.region_n] = p.Nd                  #   [Nd  0.0;
    params.doping[p.iphip, p.region_p] = p.Na                  #    0.0  Na]

    return params

end

function ChargeTransport.ParamsNodal(grid, p::tesca_setup, tesca_doping_file::String)

    paramsnodal = ChargeTransport.ParamsNodal(
        grid,
        p.numberOfCarriers
    )

    # Nodal Doping (Node dependend doping) mit Ansatz C = n-p
    TescaData = read_diodat(tesca_doping_file)

    n = TescaData["eDensity"]
    h = TescaData["hDensity"]
    NDoping = (-1) .* (n .- h) / p.cm^3 # Tesca rechnet mit cm³, ChargeTransport mit m³

    coords = grid[Coordinates]
    numberOfNodes = size(coords, 2)

    for icoord in 1:numberOfNodes
        paramsnodal.doping[icoord] = NDoping[icoord]
    end

    return paramsnodal

end
