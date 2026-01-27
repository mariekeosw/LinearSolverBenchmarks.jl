# 2D Grid
function generate_2Dgrid(grid_refinement)

    # unit factors
    @local_unitfactors μm

    # region numbers
    # Note: we need one region to create Tesca files, since the material is the same everywhere
    region_n = 1 # n-doped region
    region_p = 2 # p-doped region

    # boundary region numbers
    bregion_gate = 1
    bregion_drain = 2
    bregion_source = 3
    bregion_bulk = 4
    bregion_neumann = 5

    # Gridpoints x-direction
    x0 = -1.5 * μm              # beginning left contact (drain)
    x1 = -1.0 * μm              # end left contact (drain)
    x2 = -0.9 * μm              # beginning gate contact
    x3 = -0.8 * μm              # end n-doped zone drain
    x4 = 0.8 * μm               # beginning n-doped zone source
    x5 = 0.9 * μm               # end gate contact
    x6 = 1.0 * μm               # beginning right contact (source)
    x7 = 1.5 * μm               # end right contact (source)

    # Gridpoints y-direction
    y0 = -2.4 * μm              # bottom of the device (bulk)
    y1 = -0.5 * μm
    y2 = -0.2 * μm              # beginning n-channel (gate region)
    y3 = 0.0 * μm               # top of the device

    # Refinement in x-direction
    X = collect(range(x0, x7, length = 31 + grid_refinement * 30))

    # Refinement in y-direction
    Y1 = collect(range(y0, y1, length = 20 + grid_refinement * 19))
    Y2 = collect(range(y1, y3, length = 11 + grid_refinement * 10))
    Y = glue(Y1, Y2)

    grid = simplexgrid(X, Y)

    # cell regions
    cellmask!(grid, [x0, y2], [x3, y3], region_n)
    cellmask!(grid, [x4, y2], [x7, y3], region_n)
    cellmask!(grid, [x3, y2], [x4, y3], region_p)
    cellmask!(grid, [x0, y0], [x7, y2], region_p)

    # boundary regions
    bfacemask!(grid, [x2, y3], [x5, y3], bregion_gate)
    bfacemask!(grid, [x0, y3], [x1, y3], bregion_drain)
    bfacemask!(grid, [x6, y3], [x7, y3], bregion_source)
    bfacemask!(grid, [x0, y0], [x7, y0], bregion_bulk)
    bfacemask!(grid, [x0, y0], [x0, y3], bregion_neumann)
    bfacemask!(grid, [x7, y0], [x7, y3], bregion_neumann)
    bfacemask!(grid, [x1, y3], [x2, y3], bregion_neumann)
    bfacemask!(grid, [x5, y3], [x6, y3], bregion_neumann)

    return grid
end

# 3D Grid
function generate_3Dgrid(grid_refinement)

    # unit factors
    @local_unitfactors μm

    # region numbers
    region_n = 1 # n-doped region
    region_p = 2 # p-doped region

    # boundary region numbers
    bregion_gate = 1
    bregion_drain = 2
    bregion_source = 3
    bregion_bulk = 4
    bregion_neumann = 5

    # Gridpoints x-direction
    x0 = -1.5 * μm              # beginning left contact (drain)
    x1 = -1.0 * μm              # end left contact (drain)
    x2 = -0.9 * μm              # beginning gate contact
    x3 = -0.8 * μm              # end n-doped zone drain
    x4 = 0.8 * μm               # beginning n-doped zone source
    x5 = 0.9 * μm               # end gate contact
    x6 = 1.0 * μm               # beginning right contact (source)
    x7 = 1.5 * μm               # end right contact (source)

    # Gridpoints y-direction
    y0 = 0 * μm
    y1 = 5 * μm                # begin gate electrode
    y2 = 10 * μm               # end gate electrode
    y3 = 15 * μm               # das ist die tatsächliche Länge aus Tesca

    # Gridpoints in z-direction
    z0 = -2.4 * μm              # bottom of the device (bulk)
    z1 = -0.5 * μm
    z2 = -0.2 * μm              # beginning n-channel (gate region)
    z3 = 0.0 * μm               # top of the device

    # Refinement in x-direction
    X = collect(range(x0, x7, length = 31 + grid_refinement * 30))

    # Refinement in y-direction
    Y = collect(range(y0, y3, length = 61 + grid_refinement * 60))

    # Refinement in z-direction
    Z1 = collect(range(z0, z1, length = 20 + grid_refinement * 19))
    Z2 = collect(range(z1, z3, length = 11 + grid_refinement * 10))
    Z = glue(Z1, Z2)

    grid = simplexgrid(X, Y, Z)

    # cell regions
    cellmask!(grid, [x0, y0, z2], [x3, y3, z3], region_n)
    cellmask!(grid, [x4, y0, z2], [x7, y3, z3], region_n)
    cellmask!(grid, [x3, y0, z2], [x4, y3, z3], region_p)
    cellmask!(grid, [x0, y0, z0], [x7, y3, z2], region_p)

    # boundary regions
    bfacemask!(grid, [x2, y1, z3], [x5, y2, z3], bregion_gate)
    bfacemask!(grid, [x2, y0, z3], [x5, y1, z3], bregion_neumann)
    bfacemask!(grid, [x2, y2, z3], [x5, y3, z3], bregion_neumann)
    bfacemask!(grid, [x0, y0, z3], [x1, y3, z3], bregion_drain)
    bfacemask!(grid, [x6, y0, z3], [x7, y3, z3], bregion_source)
    bfacemask!(grid, [x0, y0, z0], [x7, y3, z0], bregion_bulk)
    bfacemask!(grid, [x0, y0, z0], [x0, y3, z3], bregion_neumann)
    bfacemask!(grid, [x7, y0, z0], [x7, y3, z3], bregion_neumann)
    bfacemask!(grid, [x1, y0, z3], [x2, y3, z3], bregion_neumann)
    bfacemask!(grid, [x5, y0, z3], [x6, y3, z3], bregion_neumann)
    bfacemask!(grid, [x0, y0, z0], [x7, y0, z3], bregion_neumann)
    bfacemask!(grid, [x0, y3, z0], [x7, y3, z3], bregion_neumann)

    return grid
end

# 3D Grid Half
function generate_3Dgrid_half(grid_refinement)

    # unit factors
    @local_unitfactors μm

    # region numbers
    region_n = 1 # n-doped region
    region_p = 2 # p-doped region

    # boundary region numbers
    bregion_gate = 1
    bregion_drain = 2
    bregion_source = 3
    bregion_bulk = 4
    bregion_neumann = 5

    # Gridpoints x-direction
    x0 = -1.5 * μm              # beginning left contact (drain)
    x1 = -1.0 * μm              # end left contact (drain)
    x2 = -0.9 * μm              # beginning gate contact
    x3 = -0.8 * μm              # end n-doped zone drain
    x4 = 0.8 * μm               # beginning n-doped zone source
    x5 = 0.9 * μm               # end gate contact
    x6 = 1.0 * μm               # beginning right contact (source)
    x7 = 1.5 * μm               # end right contact (source)

    # Gridpoints y-direction
    y0 = 0 * μm
    y1 = 5 * μm                # begin gate electrode
    y2 = 7.5 * μm              # mid point device, to use symmetry effect, take 15 * μm for full grid

    # Gridpoints in z-direction
    z0 = -2.4 * μm              # bottom of the device (bulk)
    z1 = -0.5 * μm
    z2 = -0.2 * μm              # beginning n-channel (gate region)
    z3 = 0.0 * μm               # top of the device

    # Refinement in x-direction
    X = collect(range(x0, x7, length = 31 + grid_refinement * 30))

    # Refinement in y-direction
    Y = collect(range(y0, y2, length = 31 + grid_refinement * 30))

    # Refinement in z-direction
    Z1 = collect(range(z0, z1, length = 20 + grid_refinement * 19))
    Z2 = collect(range(z1, z3, length = 11 + grid_refinement * 10))
    Z = glue(Z1, Z2)

    grid = simplexgrid(X, Y, Z)

    # cell regions
    cellmask!(grid, [x0, y0, z2], [x3, y2, z3], region_n)
    cellmask!(grid, [x4, y0, z2], [x7, y2, z3], region_n)
    cellmask!(grid, [x3, y0, z2], [x4, y2, z3], region_p)
    cellmask!(grid, [x0, y0, z0], [x7, y2, z2], region_p)

    # boundary regions
    bfacemask!(grid, [x2, y1, z3], [x5, y2, z3], bregion_gate)
    bfacemask!(grid, [x2, y0, z3], [x5, y1, z3], bregion_neumann)
    bfacemask!(grid, [x0, y0, z3], [x1, y2, z3], bregion_drain)
    bfacemask!(grid, [x6, y0, z3], [x7, y2, z3], bregion_source)
    bfacemask!(grid, [x0, y0, z0], [x7, y2, z0], bregion_bulk)
    bfacemask!(grid, [x0, y0, z0], [x0, y2, z3], bregion_neumann)
    bfacemask!(grid, [x7, y0, z0], [x7, y2, z3], bregion_neumann)
    bfacemask!(grid, [x1, y0, z3], [x2, y2, z3], bregion_neumann)
    bfacemask!(grid, [x5, y0, z3], [x6, y2, z3], bregion_neumann)
    bfacemask!(grid, [x0, y0, z0], [x7, y0, z3], bregion_neumann)
    bfacemask!(grid, [x0, y2, z0], [x7, y2, z3], bregion_neumann)

    return grid
end
