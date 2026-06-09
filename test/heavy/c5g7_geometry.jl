import Gridap: get_face_labeling
import Gridap.Geometry: get_face_tag, get_tag_from_name

const C5G7_MATERIAL_NAMES = [
    "guide-tube", "fission-chamber", "UO2", "MOX_43", "MOX_7", "MOX_87", "water"
]

const C5G7_N_PINS = 17
const C5G7_PIN_PITCH = 1.26
const C5G7_PIN_RADIUS = 0.54
const C5G7_ASSEMBLY_WIDTH = C5G7_N_PINS * C5G7_PIN_PITCH
const C5G7_CORE_WIDTH = 2 * C5G7_ASSEMBLY_WIDTH
const C5G7_REFLECTOR_FINE_PIN_LAYERS = 11
const C5G7_REFLECTOR_COARSE_PIN_LAYERS = 6
const C5G7_REFLECTOR_REFINES = 3

const C5G7_PIN_COUNT = (
    guide_tube=96,
    fission_chamber=4,
    uo2=528,
    mox_43=128,
    mox_7=200,
    mox_87=200,
)

const C5G7_ACTIVE_FUEL_PIN_COUNT =
    C5G7_PIN_COUNT.uo2 +
    C5G7_PIN_COUNT.mox_43 +
    C5G7_PIN_COUNT.mox_7 +
    C5G7_PIN_COUNT.mox_87

const C5G7_ACTIVE_FUEL_TAGS = Int32[3, 4, 5, 6]

const C5G7_ASSEMBLY_LABELS = (
    :lower_left_uo2,
    :lower_right_mox,
    :upper_left_mox,
    :upper_right_uo2,
)

function add_c5g7_circle!(factory, origin, radius, lc)
    ox, oy = origin
    center = factory.addPoint(ox, oy, 0.0, lc)
    right = factory.addPoint(ox + radius, oy, 0.0, lc)
    top = factory.addPoint(ox, oy + radius, 0.0, lc)
    left = factory.addPoint(ox - radius, oy, 0.0, lc)
    bottom = factory.addPoint(ox, oy - radius, 0.0, lc)

    c1 = factory.addCircleArc(right, center, top)
    c2 = factory.addCircleArc(top, center, left)
    c3 = factory.addCircleArc(left, center, bottom)
    c4 = factory.addCircleArc(bottom, center, right)

    return factory.addCurveLoop([c1, c2, c3, c4])
end

function c5g7_reflector_offsets()
    offsets = Float64[0.0]
    fine_width = C5G7_PIN_PITCH / C5G7_REFLECTOR_REFINES

    for i in 1:(C5G7_REFLECTOR_FINE_PIN_LAYERS * C5G7_REFLECTOR_REFINES)
        push!(offsets, i * fine_width)
    end

    fine_end = C5G7_REFLECTOR_FINE_PIN_LAYERS * C5G7_PIN_PITCH
    for i in 1:C5G7_REFLECTOR_COARSE_PIN_LAYERS
        push!(offsets, fine_end + i * C5G7_PIN_PITCH)
    end

    return offsets
end

function c5g7_refined_core_edges()
    fine_width = C5G7_PIN_PITCH / C5G7_REFLECTOR_REFINES
    n_edges = 2 * C5G7_N_PINS * C5G7_REFLECTOR_REFINES
    return [i * fine_width for i in 0:n_edges]
end

function add_c5g7_rectangular_grid!(factory, x_edges, y_edges, lc)
    points = Matrix{Int}(undef, length(x_edges), length(y_edges))
    for j in eachindex(y_edges), i in eachindex(x_edges)
        points[i, j] = factory.addPoint(x_edges[i], y_edges[j], 0.0, lc)
    end

    horizontal = Matrix{Int}(undef, length(x_edges) - 1, length(y_edges))
    for j in eachindex(y_edges), i in 1:(length(x_edges)-1)
        horizontal[i, j] = factory.addLine(points[i, j], points[i+1, j])
    end

    vertical = Matrix{Int}(undef, length(x_edges), length(y_edges) - 1)
    for j in 1:(length(y_edges)-1), i in eachindex(x_edges)
        vertical[i, j] = factory.addLine(points[i, j], points[i, j+1])
    end

    surfaces = Int[]
    sizehint!(surfaces, (length(x_edges) - 1) * (length(y_edges) - 1))
    for j in 1:(length(y_edges)-1), i in 1:(length(x_edges)-1)
        loop = factory.addCurveLoop([
            horizontal[i, j],
            vertical[i+1, j],
            -horizontal[i, j+1],
            -vertical[i, j],
        ])
        push!(surfaces, factory.addPlaneSurface([loop]))
    end

    return (;
        surfaces,
        bottom=collect(horizontal[:, 1]),
        right=collect(vertical[end, :]),
        top=collect(horizontal[:, end]),
        left=collect(vertical[1, :]),
    )
end

function c5g7_benchmark_model(; pin_lc=0.25, reflector_lc=1.0)
    dir = mktempdir()
    mshfile = joinpath(dir, "c5g7-benchmark.msh")

    gmsh.initialize()
    try
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("c5g7-benchmark")
        factory = gmsh.model.geo

        guide_tube = Int[]
        fission_chamber = Int[]
        uo2 = Int[]
        mox_43 = Int[]
        mox_7 = Int[]
        mox_87 = Int[]
        pin_loops = Int[]
        n_pins = C5G7_N_PINS
        pitch = C5G7_PIN_PITCH
        radius = C5G7_PIN_RADIUS

        gt_pos = [(3, 6), (3, 9), (3, 12),
                  (4, 4), (4, 14),
                  (6, 3), (6, 6), (6, 9), (6, 12), (6, 15),
                  (9, 3), (9, 6), (9, 12), (9, 15),
                  (12, 3), (12, 6), (12, 9), (12, 12), (12, 15),
                  (14, 4), (14, 14),
                  (15, 6), (15, 9), (15, 12)]
        fc_pos = [(9, 9)]

        for n in 1:2
            x0 = (n - 1) * n_pins * pitch
            y0 = (n - 1) * n_pins * pitch
            for i in 1:n_pins, j in 1:n_pins
                x = x0 + pitch / 2 + (i - 1) * pitch
                y = y0 + pitch / 2 + (j - 1) * pitch
                loop = add_c5g7_circle!(factory, (x, y), radius, pin_lc)
                surface = factory.addPlaneSurface([loop])
                push!(pin_loops, loop)

                pos = (i, j)
                if pos in gt_pos
                    push!(guide_tube, surface)
                elseif pos in fc_pos
                    push!(fission_chamber, surface)
                else
                    push!(uo2, surface)
                end
            end
        end

        mox_43_pos = Tuple{Int,Int}[]
        for i in 1:n_pins, j in 1:n_pins
            if i == 1 || i == n_pins || j == 1 || j == n_pins
                push!(mox_43_pos, (i, j))
            end
        end

        mox_87_pos = Tuple{Int,Int}[]
        for i in 4:(n_pins-3), j in 4:(n_pins-3)
            pos = (i, j)
            if pos in gt_pos || pos in fc_pos ||
               pos in ((5, 4), (13, 4), (4, 5), (14, 5),
                       (4, 13), (14, 13), (5, 14), (13, 14))
                continue
            end
            push!(mox_87_pos, pos)
        end

        mox_7_pos = Tuple{Int,Int}[]
        for i in 1:n_pins, j in 1:n_pins
            pos = (i, j)
            if pos in gt_pos || pos in fc_pos || pos in mox_43_pos || pos in mox_87_pos
                continue
            end
            push!(mox_7_pos, pos)
        end

        for n in 1:2
            x0 = isone(n) ? n_pins * pitch : 0.0
            y0 = isone(n) ? 0.0 : n_pins * pitch
            for i in 1:n_pins, j in 1:n_pins
                x = x0 + pitch / 2 + (i - 1) * pitch
                y = y0 + pitch / 2 + (j - 1) * pitch
                loop = add_c5g7_circle!(factory, (x, y), radius, pin_lc)
                surface = factory.addPlaneSurface([loop])
                push!(pin_loops, loop)

                pos = (i, j)
                if pos in gt_pos
                    push!(guide_tube, surface)
                elseif pos in mox_43_pos
                    push!(mox_43, surface)
                elseif pos in mox_7_pos
                    push!(mox_7, surface)
                elseif pos in mox_87_pos
                    push!(mox_87, surface)
                elseif pos in fc_pos
                    push!(fission_chamber, surface)
                end
            end
        end

        core_side = 2 * n_pins * pitch
        p1 = factory.addPoint(0.0, 0.0, 0.0, pin_lc)
        p2 = factory.addPoint(core_side, 0.0, 0.0, pin_lc)
        p3 = factory.addPoint(core_side, core_side, 0.0, pin_lc)
        p4 = factory.addPoint(0.0, core_side, 0.0, pin_lc)

        core_bottom = factory.addLine(p1, p2)
        core_right = factory.addLine(p2, p3)
        core_top = factory.addLine(p3, p4)
        core_left = factory.addLine(p4, p1)
        core_loop = factory.addCurveLoop([core_bottom, core_right, core_top, core_left])
        moderator = factory.addPlaneSurface(
            [core_loop; pin_loops]
        )

        core_edges = c5g7_refined_core_edges()
        reflector_offsets = c5g7_reflector_offsets()
        reflector_edges = core_side .+ reflector_offsets

        right_reflector = add_c5g7_rectangular_grid!(
            factory, reflector_edges, core_edges, reflector_lc
        )
        top_reflector = add_c5g7_rectangular_grid!(
            factory, core_edges, reflector_edges, reflector_lc
        )
        corner_reflector = add_c5g7_rectangular_grid!(
            factory, reflector_edges, reflector_edges, reflector_lc
        )

        reflector_surfaces = [
            right_reflector.surfaces;
            top_reflector.surfaces;
            corner_reflector.surfaces
        ]

        surface_groups = (
            guide_tube, fission_chamber, uo2, mox_43, mox_7, mox_87,
            [moderator; reflector_surfaces]
        )
        for (name, surfaces) in zip(C5G7_MATERIAL_NAMES, surface_groups)
            group = factory.addPhysicalGroup(2, surfaces)
            gmsh.model.setPhysicalName(2, group, name)
        end

        boundary_groups = (
            bottom=[core_bottom],
            left=[core_left],
            top=[top_reflector.top; corner_reflector.top],
            right=[right_reflector.right; corner_reflector.right],
        )
        for (name, lines) in pairs(boundary_groups)
            group = factory.addPhysicalGroup(1, lines)
            gmsh.model.setPhysicalName(1, group, String(name))
        end

        factory.removeAllDuplicates()
        factory.synchronize()
        gmsh.model.mesh.generate(2)
        gmsh.write(mshfile)
    finally
        gmsh.finalize()
    end

    return GmshDiscreteModel(mshfile; renumber=true)
end

function c5g7_analytic_material_volumes()
    pin_area = pi * C5G7_PIN_RADIUS^2
    total_area = (3 * C5G7_ASSEMBLY_WIDTH)^2

    guide_tube = C5G7_PIN_COUNT.guide_tube * pin_area
    fission_chamber = C5G7_PIN_COUNT.fission_chamber * pin_area
    uo2 = C5G7_PIN_COUNT.uo2 * pin_area
    mox_43 = C5G7_PIN_COUNT.mox_43 * pin_area
    mox_7 = C5G7_PIN_COUNT.mox_7 * pin_area
    mox_87 = C5G7_PIN_COUNT.mox_87 * pin_area
    water = total_area - guide_tube - fission_chamber - uo2 - mox_43 - mox_7 - mox_87

    return [guide_tube, fission_chamber, uo2, mox_43, mox_7, mox_87, water]
end

function c5g7_material_volumes(model)
    mesh = RayTracing.Mesh(model)
    face_labeling = get_face_labeling(model)
    gridap_tags = [get_tag_from_name(face_labeling, name) for name in C5G7_MATERIAL_NAMES]
    tag_to_idx = Dict(Int(tag) => i for (i, tag) in enumerate(gridap_tags))
    cell_tags = get_face_tag(face_labeling, C5G7_MATERIAL_NAMES, 2)

    volumes = zeros(length(C5G7_MATERIAL_NAMES))
    for (cell, tag) in enumerate(cell_tags)
        volumes[tag_to_idx[Int(tag)]] += RayTracing.element_volume(mesh, cell)
    end

    return volumes
end

function c5g7_cell_centroid(mesh, cell::Integer)
    node_ids = mesh.ordered_cell_nodes[cell]
    x = 0.0
    y = 0.0
    for node_id in node_ids
        node = mesh.node_coordinates[node_id]
        x += node[1]
        y += node[2]
    end
    n = length(node_ids)
    return x / n, y / n
end

function c5g7_assembly_index(x::Real, y::Real)
    if !(0 <= x <= C5G7_CORE_WIDTH && 0 <= y <= C5G7_CORE_WIDTH)
        return nothing
    end

    col = clamp(floor(Int, x / C5G7_ASSEMBLY_WIDTH) + 1, 1, 2)
    row = clamp(floor(Int, y / C5G7_ASSEMBLY_WIDTH) + 1, 1, 2)

    if col == 1 && row == 1
        return 1
    elseif col == 2 && row == 1
        return 2
    elseif col == 1 && row == 2
        return 3
    elseif col == 2 && row == 2
        return 4
    end

    return nothing
end

function c5g7_pin_index(x::Real, y::Real)
    assembly_idx = c5g7_assembly_index(x, y)
    isnothing(assembly_idx) && return nothing

    col = isodd(assembly_idx) ? 1 : 2
    row = assembly_idx <= 2 ? 1 : 2
    local_x = x - (col - 1) * C5G7_ASSEMBLY_WIDTH
    local_y = y - (row - 1) * C5G7_ASSEMBLY_WIDTH

    pin_i = clamp(floor(Int, local_x / C5G7_PIN_PITCH) + 1, 1, C5G7_N_PINS)
    pin_j = clamp(floor(Int, local_y / C5G7_PIN_PITCH) + 1, 1, C5G7_N_PINS)

    return assembly_idx, pin_i, pin_j
end
