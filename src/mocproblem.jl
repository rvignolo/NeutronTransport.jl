import RayTracing: num_dims, num_cells

"""
    MoCProblem{Dim,NRegions,NGroups,G<:TrackGenerator,Q<:Quadrature,M} <: TransportProblem{Dim,NRegions,NGroups}

Defines a problem that ought to be solved using the Method of Characteristics, with
dimension `Dim`, number of flat source regions `NRegions`, number of energy groups
`NGroups`, ray tracing
"""
struct MoCProblem{Dim,NRegions,NGroups,G<:TrackGenerator,Q<:Quadrature,Xs<:XSs} <: TransportProblem{Dim,NRegions,NGroups}
    nφ::Int
    nψ::Int

    trackgenerator::G
    quadrature::Q

    xss::Vector{Xs}
    fsr_tag::Vector{Int8}  # cell/element/fsr to xs tag (tag is the index in the xss array)

    function MoCProblem{Dim,NRegions,NGroups}(
        nφ, nψ, tg::G, quad::Q, xss::Vector{X}, fsr_tag
    ) where {Dim,NRegions,NGroups,G,Q,X}
        return new{Dim,NRegions,NGroups,G,Q,X}(nφ, nψ, tg, quad, xss, fsr_tag)
    end
end

function MoCProblem(tg::TrackGenerator, polar_quadrature::PolarQuadrature, xss::Vector{<:XSs})
    @unpack mesh, azimuthal_quadrature, n_total_tracks = tg

    Dim = num_dims(mesh)

    NGroups = ngroups(first(xss))
    if !all(g -> isequal(g, NGroups), ngroups.(xss))
        error("all `CrossSections` *must* have the same number of energy groups.")
    end

    NRegions = num_cells(mesh)
    n_polar_2 = npolar2(polar_quadrature)

    nφ = NRegions * NGroups
    nψ = n_total_tracks * 2 * n_polar_2 * NGroups # 2 is the number of directions

    quadrature = Quadrature(azimuthal_quadrature, polar_quadrature)

    face_labeling = get_face_labeling(tg.mesh.model)
    tag_to_idx = Dict{Int8,Int8}()
    for (i, xs) in enumerate(xss)
        # given a xs name, find its gridap tag
        gridap_tag = get_tag_from_name(face_labeling, xs.name)

        # map gridap tag to xs idx
        push!(tag_to_idx, gridap_tag => i)
    end

    cell_tag = get_face_tag(face_labeling, Dim)  # gridap cell/element/fsr tags
    fsr_tag = similar(cell_tag)                  # NeutronTransport cell/element/fsr tags
    for (i, gridap_tag) in enumerate(cell_tag)
        fsr_tag[i] = tag_to_idx[gridap_tag]
    end

    return MoCProblem{Dim,NRegions,NGroups}(nφ, nψ, tg, quadrature, xss, fsr_tag)
end

dimension(::MoCProblem{Dim}) where{Dim} = Dim
nregions(::MoCProblem{Dim,NRegions}) where{Dim,NRegions} = NRegions
ngroups(::MoCProblem{Dim,NRegions,NGroups}) where{Dim,NRegions,NGroups} = NGroups

function show(io::IO, prob::MoCProblem)
    println(io, "  Problem Dimension: ", dimension(prob))
    println(io, "  Number of regions: ", nregions(prob))
    print(io,   "  Number of energy groups: ", ngroups(prob))
end

# si quiero que los primeros NRegions sean para el 1er group, etc...
# macro region_index(i, g)
#     return :((g - 1) * NRegions + i)
# end

# si quiero que esten intercalados
macro region_index(i, g)
    ex = quote
        ($i - 1) * NGroups + $g
    end
    return esc(ex)
end
# puedo hacerlas como function con @inline pero necesitan `groups` o `nfsr`

macro angular_index(t, d, p, g)
    ex = quote
        ($t - 1) * 2 * n_polar_2 * NGroups + $d * n_polar_2 * NGroups + ($p - 1) * NGroups + $g
    end
    return esc(ex)
end

macro reduced_angular_index(p, g)
    ex = quote
        ($p - 1) * NGroups + $g
    end
    return esc(ex)
end
