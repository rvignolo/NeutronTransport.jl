
struct NeutronTransportProblem{F,S}
    formulation::F
    solver::S
end

NeutronTransportProblem(
    f::MethodOfCharacteristics, tg::TrackGenerator, pq::PolarQuadrature
) = NeutronTransportProblem(f, MocSolver(tg, pq))

NeutronTransportProblem(f::CollisionProbability) = error("not yet implemented.")