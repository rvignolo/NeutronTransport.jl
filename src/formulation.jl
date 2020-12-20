abstract type Formulation end

struct MethodOfCharacteristics <: Formulation end
struct CollisionProbability <: Formulation end
struct DiscreteOrdinates <: Formulation end
struct Diffusion <: Formulation end

const MoC = MethodOfCharacteristics
const CP = CollisionProbability
const SN = DiscreteOrdinates
