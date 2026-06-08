import Gridap: get_face_labeling
import Gridap.Geometry: get_face_tag, get_tag_from_name

const C5G7_MATERIAL_NAMES = [
    "guide-tube", "fission-chamber", "UO2", "MOX_43", "MOX_7", "MOX_87", "water"
]

const C5G7_PIN_COUNT = (
    guide_tube=96,
    fission_chamber=4,
    uo2=528,
    mox_43=128,
    mox_7=200,
    mox_87=200,
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

function c5g7_benchmark_model(; pin_lc=0.25, reflector_lcs=(1.0, 1.5, 2.0))
    dir = mktempdir()
    mshfile = joinpath(dir, "c5g7-benchmark.msh")

    n_pins = 17
    pitch = 1.26
    radius = 0.54

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

        reflector_surfaces = Int[]
        bottom_lines = [core_bottom]
        left_lines = [core_left]
        top_line = core_top
        right_line = core_right
        side = core_side
        reflector_width = n_pins * pitch / 3

        for lc in reflector_lcs
            p1 = factory.addPoint(side, 0.0, 0.0, lc)
            p2 = factory.addPoint(side, side, 0.0, lc)
            p3 = factory.addPoint(0.0, side, 0.0, lc)
            p4 = factory.addPoint(0.0, side + reflector_width, 0.0, lc)
            p5 = factory.addPoint(side + reflector_width, side + reflector_width, 0.0, lc)
            p6 = factory.addPoint(side + reflector_width, 0.0, 0.0, lc)

            inner_right = factory.addLine(p1, p2)
            inner_top = factory.addLine(p2, p3)
            left_extension = factory.addLine(p3, p4)
            outer_top = factory.addLine(p4, p5)
            outer_right = factory.addLine(p5, p6)
            bottom_extension = factory.addLine(p6, p1)

            loop = factory.addCurveLoop([
                inner_right, inner_top, left_extension,
                outer_top, outer_right, bottom_extension
            ])
            push!(reflector_surfaces, factory.addPlaneSurface([loop]))
            push!(left_lines, left_extension)
            push!(bottom_lines, bottom_extension)
            top_line = outer_top
            right_line = outer_right
            side += reflector_width
        end

        surface_groups = (
            guide_tube, fission_chamber, uo2, mox_43, mox_7, mox_87,
            [moderator; reflector_surfaces]
        )
        for (name, surfaces) in zip(C5G7_MATERIAL_NAMES, surface_groups)
            group = factory.addPhysicalGroup(2, surfaces)
            gmsh.model.setPhysicalName(2, group, name)
        end

        boundary_groups = (
            bottom=bottom_lines,
            left=left_lines,
            top=[top_line],
            right=[right_line],
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
    n_pins = 17
    pitch = 1.26
    radius = 0.54
    pin_area = pi * radius^2
    total_area = (3 * n_pins * pitch)^2

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
