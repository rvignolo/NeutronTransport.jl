import Gridap: DiscreteModelFromFile
import Gridap.Geometry:
    CartesianDiscreteModel,
    FaceLabeling,
    UnstructuredDiscreteModel,
    get_grid,
    get_grid_topology
import GridapGmsh: GmshDiscreteModel, gmsh

const TEST_REPO_ROOT = dirname(@__DIR__)
const TEST_DEMO_DIR = joinpath(TEST_REPO_ROOT, "demo")

function single_cell_model()
    cartesian = CartesianDiscreteModel((0.0, 1.0, 0.0, 1.0), (1, 1))
    base = UnstructuredDiscreteModel(cartesian)
    topology = get_grid_topology(base)
    labels = FaceLabeling(topology, [1], ["fuel"])
    return UnstructuredDiscreteModel(get_grid(base), topology, labels)
end

function cartesian_labeled_model(domain, partition, cell_to_tag, tag_to_name)
    cartesian = CartesianDiscreteModel(domain, partition)
    base = UnstructuredDiscreteModel(cartesian)
    topology = get_grid_topology(base)
    labels = FaceLabeling(topology, cell_to_tag, tag_to_name)
    return UnstructuredDiscreteModel(get_grid(base), topology, labels)
end

function single_cell_problem(bcs)
    tg = TrackGenerator(single_cell_model(), 4, 0.3; bcs=bcs)
    trace!(tg)
    segmentize!(tg)

    fuel = CrossSections("fuel", 2;
        νΣf=[0.02, 0.3],
        Σt=[0.4, 0.7],
        Σs0=[0.3 0.05; 0.0 0.4]
    )

    prob = MoCProblem(tg, TabuchiYamamoto(2), [fuel])
    return prob, tg
end

function solve_single_cell(bcs)
    prob, tg = single_cell_problem(bcs)
    sol = NeutronTransport.solve(prob; max_iterations=3, max_residual=0.0)

    return prob, sol, tg
end

function openmoc_homogeneous_medium()
    return CrossSections("fuel", 2;
        νΣf=[0.0015, 0.325],
        χ=[1.0, 0.0],
        Σt=[0.2208, 1.604],
        Σs0=[0.1 0.117; 0.0 1.42]
    )
end

function openmoc_homogeneous_medium_keff()
    xs = openmoc_homogeneous_medium()
    flux_ratio = xs.Σs0[1, 2] / (xs.Σt[2] - xs.Σs0[2, 2])
    keff = (xs.νΣf[1] + xs.νΣf[2] * flux_ratio) / (xs.Σt[1] - xs.Σs0[1, 1])
    return keff, flux_ratio
end

function openmoc_homogeneous_grid_model()
    length = 2.5
    num_cells = 10
    cell_to_tag = fill(1, num_cells * num_cells)
    return cartesian_labeled_model(
        (-length / 2, length / 2, -length / 2, length / 2),
        (num_cells, num_cells),
        cell_to_tag,
        ["fuel"]
    )
end

function openmoc_reflective_grid_model()
    num_cells = 3
    cell_to_tag = fill(2, num_cells * num_cells)
    center_cell = div(length(cell_to_tag), 2) + 1
    cell_to_tag[center_cell] = 1
    return cartesian_labeled_model(
        (-3.0, 3.0, -3.0, 3.0),
        (num_cells, num_cells),
        cell_to_tag,
        ["fuel", "moderator"]
    )
end

function openmoc_lattice_grid_model()
    # OpenMOC represents this with a 3x3 lattice. In Gridap we use the equivalent
    # labeled Cartesian mesh: one central fuel FSR surrounded by moderator FSRs.
    cell_to_tag = [
        2, 2, 2,
        2, 1, 2,
        2, 2, 2,
    ]
    return cartesian_labeled_model(
        (-3.0, 3.0, -3.0, 3.0),
        (3, 3),
        cell_to_tag,
        ["fuel", "moderator"]
    )
end

function openmoc_oblong_lattice_grid_model()
    # OpenMOC's oblong lattice has this row-major material pattern:
    # [fuel guide-tube moderator; moderator fuel guide-tube].
    cell_to_tag = [1, 3, 2, 2, 1, 3]
    return cartesian_labeled_model(
        (-3.0, 3.0, -2.0, 2.0),
        (3, 2),
        cell_to_tag,
        ["fuel", "moderator", "guide-tube"]
    )
end

function openmoc_water_medium()
    return CrossSections("moderator", 2;
        Σt=[0.640711, 1.69131],
        Σs0=[0.607382 0.0331316; 0.0 1.68428]
    )
end

function openmoc_guide_tube_medium()
    return CrossSections("guide-tube", 2;
        Σt=[0.274144, 0.280890],
        Σs0=[0.272377 0.000190838; 0.0 0.277230]
    )
end

function openmoc_pin_cell_model(; lc=0.18)
    dir = mktempdir()
    mshfile = joinpath(dir, "openmoc-pin-cell.msh")

    gmsh.initialize()
    try
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("openmoc-pin-cell")
        factory = gmsh.model.geo

        # OpenMOC PinCellInput: one fuel cylinder of radius 1.0 centered at the
        # origin, bounded by a reflective square from -2 to 2 in x and y.
        center = factory.addPoint(0.0, 0.0, 0.0, lc)
        right = factory.addPoint(1.0, 0.0, 0.0, lc)
        top = factory.addPoint(0.0, 1.0, 0.0, lc)
        left = factory.addPoint(-1.0, 0.0, 0.0, lc)
        bottom = factory.addPoint(0.0, -1.0, 0.0, lc)

        c1 = factory.addCircleArc(right, center, top)
        c2 = factory.addCircleArc(top, center, left)
        c3 = factory.addCircleArc(left, center, bottom)
        c4 = factory.addCircleArc(bottom, center, right)
        fuel_loop = factory.addCurveLoop([c1, c2, c3, c4])

        p1 = factory.addPoint(-2.0, -2.0, 0.0, lc)
        p2 = factory.addPoint(2.0, -2.0, 0.0, lc)
        p3 = factory.addPoint(2.0, 2.0, 0.0, lc)
        p4 = factory.addPoint(-2.0, 2.0, 0.0, lc)

        bottom_line = factory.addLine(p1, p2)
        right_line = factory.addLine(p2, p3)
        top_line = factory.addLine(p3, p4)
        left_line = factory.addLine(p4, p1)
        moderator_loop = factory.addCurveLoop([bottom_line, right_line, top_line, left_line])

        fuel = factory.addPlaneSurface([fuel_loop])
        moderator = factory.addPlaneSurface([moderator_loop, fuel_loop])

        fuel_group = factory.addPhysicalGroup(2, [fuel])
        moderator_group = factory.addPhysicalGroup(2, [moderator])
        bottom_group = factory.addPhysicalGroup(1, [bottom_line])
        right_group = factory.addPhysicalGroup(1, [right_line])
        top_group = factory.addPhysicalGroup(1, [top_line])
        left_group = factory.addPhysicalGroup(1, [left_line])

        gmsh.model.setPhysicalName(2, fuel_group, "UO2")
        gmsh.model.setPhysicalName(2, moderator_group, "Water")
        gmsh.model.setPhysicalName(1, bottom_group, "bottom")
        gmsh.model.setPhysicalName(1, right_group, "right")
        gmsh.model.setPhysicalName(1, top_group, "top")
        gmsh.model.setPhysicalName(1, left_group, "left")

        factory.synchronize()
        gmsh.model.mesh.generate(2)
        gmsh.write(mshfile)
    finally
        gmsh.finalize()
    end

    return GmshDiscreteModel(mshfile; renumber=true)
end

function openmoc_pin_cell_materials()
    n_groups = 7
    χ = round.([0.58791, 0.41176, 0.00033906, 1.1761e-7, 0.0, 0.0, 0.0]; digits=4)

    uo2 = CrossSections("UO2", n_groups;
        χ=χ,
        νΣf=[0.02005998, 0.002027303, 0.01570599, 0.04518301, 0.04334208, 0.2020901, 0.5257105],
        Σt=[0.177949, 0.329805, 0.480388, 0.554367, 0.311801, 0.395168, 0.564406],
        Σs0=[0.127537 0.042378 9.4374e-6 5.5163e-9 0.0 0.0 0.0;
            0.0 0.324456 0.0016314 3.1427e-9 0.0 0.0 0.0;
            0.0 0.0 0.45094 0.0026792 0.0 0.0 0.0;
            0.0 0.0 0.0 0.452565 0.0055664 0.0 0.0;
            0.0 0.0 0.0 0.00012525 0.271401 0.010255 1.0021e-8;
            0.0 0.0 0.0 0.0 0.0012968 0.265802 0.016809;
            0.0 0.0 0.0 0.0 0.0 0.0085458 0.27308]
    )

    water = CrossSections("Water", n_groups;
        Σt=[0.159206, 0.41297, 0.59031, 0.58435, 0.718, 1.25445, 2.65038],
        Σs0=[0.0444777 0.1134 0.00072347 3.7499e-6 5.3184e-8 0.0 0.0;
            0.0 0.282334 0.12994 0.0006234 4.8002e-5 7.4486e-6 1.0455e-6;
            0.0 0.0 0.345256 0.22457 0.016999 0.0026443 0.00050344;
            0.0 0.0 0.0 0.0910284 0.41551 0.063732 0.012139;
            0.0 0.0 0.0 7.1437e-5 0.139138 0.51182 0.061229;
            0.0 0.0 0.0 0.0 0.0022157 0.699913 0.53732;
            0.0 0.0 0.0 0.0 0.0 0.13244 2.4807]
    )

    return [uo2, water]
end

function openmoc_pin_cell_problem(;
    lc=0.18, n_azim=4, spacing=0.1, n_polar=6, coalesce_materials=false
)
    bcs = BoundaryConditions(
        top=Reflective, bottom=Reflective, left=Reflective, right=Reflective
    )
    tg = TrackGenerator(
        openmoc_pin_cell_model(; lc), n_azim, spacing; bcs=bcs, volume_correction=true
    )
    trace!(tg)
    segmentize!(tg)

    materials = openmoc_pin_cell_materials()
    cell_to_fsr = coalesce_materials ? cell_material_ids(tg, materials) : nothing
    return MoCProblem(tg, TabuchiYamamoto(n_polar), materials; cell_to_fsr=cell_to_fsr), tg
end

function demo_model(name)
    return DiscreteModelFromFile(joinpath(TEST_DEMO_DIR, "$name.json"))
end

function pincell_materials()
    pin = CrossSections("pin", 2;
        νΣf=[1.86278e-2, 3.44137e-1],
        Σt=[3.62022e-1, 5.72155e-1],
        Σs0=[3.33748e-1 6.64881e-4; 0.0e-0 3.80898e-1]
    )

    cladding = CrossSections("cladding", 2;
        Σt=[2.74144e-1, 2.80890e-1],
        Σs0=[2.72377e-1 1.90838e-4; 0.0e-0 2.77230e-1]
    )

    water = CrossSections("water", 2;
        Σt=[6.40711e-1, 1.69131e-0],
        Σs0=[6.07382e-1 3.31316e-2; 0.0e-0 1.68428e-0]
    )

    return [pin, cladding, water]
end

function bwr_materials()
    pin = CrossSections("pin", 2;
        νΣf=[1.86278e-2, 3.44137e-1],
        Σt=[3.62022e-1, 5.72155e-1],
        Σs0=[3.33748e-1 6.64881e-4; 0.0e-0 3.80898e-1]
    )

    cladding = CrossSections("cladding", 2;
        Σt=[2.74144e-1, 2.80890e-1],
        Σs0=[2.72377e-1 1.90838e-4; 0.0e-0 2.77230e-1]
    )

    water = CrossSections("water", 2;
        Σt=[6.40711e-1, 1.69131e-0],
        Σs0=[6.07382e-1 3.31316e-2; 0.0e-0 1.68428e-0]
    )

    pin_gd = CrossSections("pin-gd", 2;
        νΣf=[1.79336e-2, 1.57929e-1],
        Σt=[3.71785e-1, 1.75000e-0],
        Σs0=[3.38096e-1 6.92807e-4; 0.0e-0 3.83204e-1]
    )

    return [pin, cladding, water, pin_gd]
end

function c5g7_geometry_test_materials()
    guide_tube = CrossSections("guide-tube", 2;
        Σt=[0.274144, 0.280890],
        Σs0=[0.272377 0.000190838; 0.0 0.277230]
    )

    fission_chamber = CrossSections("fission-chamber", 2;
        Σt=[0.274144, 0.280890],
        Σs0=[0.272377 0.000190838; 0.0 0.277230]
    )

    uo2 = CrossSections("UO2", 2;
        νΣf=[1.86278e-2, 3.44137e-1],
        Σt=[3.62022e-1, 5.72155e-1],
        Σs0=[3.33748e-1 6.64881e-4; 0.0e-0 3.80898e-1]
    )

    mox_43 = CrossSections("MOX_43", 2;
        νΣf=[2.17530e-2, 6.66651e-1],
        Σt=[3.78731e-1, 6.78997e-1],
        Σs0=[3.28876e-1 8.229e-4; 0.0e-0 4.26997e-1]
    )

    mox_7 = CrossSections("MOX_7", 2;
        νΣf=[2.381395e-2, 9.281814e-1],
        Σt=[3.81323e-1, 8.33601e-1],
        Σs0=[3.30457e-1 8.5105e-4; 0.0e-0 4.74198e-1]
    )

    mox_87 = CrossSections("MOX_87", 2;
        νΣf=[2.5186e-2, 1.074999],
        Σt=[3.83045e-1, 9.21028e-1],
        Σs0=[3.31504e-1 8.6972e-4; 0.0e-0 5.02754e-1]
    )

    water = CrossSections("water", 2;
        Σt=[6.40711e-1, 1.69131e-0],
        Σs0=[6.07382e-1 3.31316e-2; 0.0e-0 1.68428e-0]
    )

    return [guide_tube, fission_chamber, uo2, mox_43, mox_7, mox_87, water]
end

function reflected_problem(model, materials; n_azim=8, spacing=0.25, n_polar=2)
    bcs = BoundaryConditions(
        top=Reflective, bottom=Reflective, left=Reflective, right=Reflective
    )
    tg = TrackGenerator(model, n_azim, spacing; bcs=bcs, volume_correction=true)
    trace!(tg)
    segmentize!(tg)
    return MoCProblem(tg, TabuchiYamamoto(n_polar), materials)
end

function demo_problem(name, materials; n_azim, spacing, n_polar=2, bcs)
    tg = TrackGenerator(demo_model(name), n_azim, spacing; bcs=bcs, volume_correction=true)
    trace!(tg)
    segmentize!(tg)
    return MoCProblem(tg, TabuchiYamamoto(n_polar), materials), tg
end
