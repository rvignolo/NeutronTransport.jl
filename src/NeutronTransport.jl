module NeutronTransport

using UnPack
using Setfield
using Reexport
using StaticArrays
@reexport using RayTracing
using Gridap: get_face_labeling
using Gridap.Geometry: get_face_tag, get_tag_from_name

import Base: show

include("formulation.jl")
include("materials.jl")
include("polar_quad.jl")
include("quadrature.jl")

abstract type TransportProblem{Dim,NRegions,NGroups} end

dimension(::TransportProblem{Dim}) where {Dim} = Dim
nregions(::TransportProblem{Dim,NRegions}) where {Dim,NRegions} = NRegions
ngroups(::TransportProblem{Dim,NRegions,NGroups}) where {Dim,NRegions,NGroups} = NGroups

include("mocproblem.jl")

abstract type TransportSolution end

include("mocsolver.jl")

export TabuchiYamamoto, GaussLegendre
export CrossSections
export MoCProblem
export solve

end