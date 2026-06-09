using NeutronTransport

import Gridap: DiscreteModelFromFile
import Gridap.Geometry:
    CartesianDiscreteModel,
    FaceLabeling,
    UnstructuredDiscreteModel,
    get_grid,
    get_grid_topology

const REPO_ROOT = dirname(@__DIR__)
const DEMO_DIR = joinpath(REPO_ROOT, "demo")

function benchmark_single_cell_model()
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

function benchmark_single_cell_material()
    return CrossSections("fuel", 2;
        νΣf=[0.02, 0.3],
        Σt=[0.4, 0.7],
        Σs0=[0.3 0.05; 0.0 0.4]
    )
end

function openmoc_homogeneous_medium()
    return CrossSections("fuel", 2;
        νΣf=[0.0015, 0.325],
        χ=[1.0, 0.0],
        Σt=[0.2208, 1.604],
        Σs0=[0.1 0.117; 0.0 1.42]
    )
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

function openmoc_water_medium()
    return CrossSections("moderator", 2;
        Σt=[0.640711, 1.69131],
        Σs0=[0.607382 0.0331316; 0.0 1.68428]
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
    cell_to_tag = [1, 3, 2, 2, 1, 3]
    return cartesian_labeled_model(
        (-3.0, 3.0, -2.0, 2.0),
        (3, 2),
        cell_to_tag,
        ["fuel", "moderator", "guide-tube"]
    )
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

function openmoc_guide_tube_medium()
    return CrossSections("guide-tube", 2;
        Σt=[0.274144, 0.280890],
        Σs0=[0.272377 0.000190838; 0.0 0.277230]
    )
end

function benchmark_single_cell_problem(;
    bcs=BoundaryConditions(),
    material=benchmark_single_cell_material(),
    n_azim=4,
    spacing=0.3,
    n_polar=2,
)
    tg = TrackGenerator(benchmark_single_cell_model(), n_azim, spacing; bcs=bcs)
    trace!(tg)
    segmentize!(tg)
    prob = MoCProblem(tg, TabuchiYamamoto(n_polar), [material])
    return prob
end

function initialized_solution(prob)
    sol = NeutronTransport.MoCSolution{eltype(prob)}(prob)
    NeutronTransport.optical_length!(prob)
    NeutronTransport.set_uniform_φ!(sol, one(eltype(prob)))
    NeutronTransport.set_uniform_Q!(sol, zero(eltype(prob)))
    NeutronTransport.set_uniform_start_boundary_ψ!(sol, zero(eltype(prob)))
    NeutronTransport.update_boundary_ψ!(sol)
    NeutronTransport.normalize_fluxes!(sol, prob)
    NeutronTransport.compute_q!(sol, prob)
    return sol
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

function demo_model(name)
    return DiscreteModelFromFile(joinpath(DEMO_DIR, "$name.json"))
end

function pincell_problem(; n_azim=8, spacing=0.05, n_polar=2)
    bcs = BoundaryConditions(
        top=Reflective, bottom=Reflective, left=Reflective, right=Reflective
    )
    tg = TrackGenerator(demo_model("pincell"), n_azim, spacing; bcs=bcs, volume_correction=true)
    trace!(tg)
    segmentize!(tg)
    return MoCProblem(tg, TabuchiYamamoto(n_polar), pincell_materials())
end

function bwr_problem(; n_azim=4, spacing=0.18, n_polar=2)
    bcs = BoundaryConditions(
        top=Reflective, bottom=Reflective, left=Reflective, right=Reflective
    )
    tg = TrackGenerator(demo_model("bwr"), n_azim, spacing; bcs=bcs, volume_correction=true)
    trace!(tg)
    segmentize!(tg)
    return MoCProblem(tg, TabuchiYamamoto(n_polar), bwr_materials())
end

function c5g7_geometry_problem(; n_azim=4, spacing=0.5, n_polar=2)
    bcs = BoundaryConditions(top=Vacuum, bottom=Reflective, left=Reflective, right=Vacuum)
    tg = TrackGenerator(demo_model("c5g7"), n_azim, spacing; bcs=bcs, volume_correction=true)
    trace!(tg)
    segmentize!(tg)
    return MoCProblem(tg, TabuchiYamamoto(n_polar), c5g7_geometry_test_materials())
end
