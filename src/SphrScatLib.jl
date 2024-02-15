module SphrScatLib
import Pkg
using LinearAlgebra
include("compute_η.jl")
include("compute_U.jl")
include("compute_coefficient_lower_hemisphere.jl")
include("compute_coefficient_out_of_hemisphere.jl")
include("Sphere.jl")
include("compute_coefficient_upper_hemishpere.jl")
include("compute_E_PEC.jl")
include("Sphere_with_Core.jl")
include("Double_hemisphere.jl")

export calculate_η
export generate_U_matrix
export calculate_LH_A2
export calculate_OH_E
export calculate_in
export calculate_HH_A1
export compute_E_PEC
export calculate_BSphere
export calculate_ASph_with_Core
export calculate_ESph_with_Core
export calculate_B
export calculate_C
export calculate_D
end # module SphrScatLib


