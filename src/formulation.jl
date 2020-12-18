abstract type Formulation end

struct MethodOfCharacteristics <: Formulation end
struct CollisionProbability <: Formulation end

const MoC = MethodOfCharacteristics
const CP = CollisionProbability