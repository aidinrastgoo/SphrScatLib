module SphrScatLib
import Pkg
using LinearAlgebra
using IterativeSolvers
using SpecialFunctions  # For Legendre polynomials
using LegendrePolynomials

include("compute_η.jl")
include("compute_U.jl")
include("compute_coefficient_lower_hemisphere.jl")
include("compute_Potential.jl")
include("Sphere.jl")
include("compute_coefficient_upper_hemishpere.jl")
include("Sphere_with_Core.jl")
include("Double_hemisphere.jl")

export calculate_η
export generate_U_matrix
export calculate_LH_A2
export calculate_in
export calculate_HH_A1
export coefficients
export calculate_Phi
export calculate_Phi_1
export calculate_Phi_2
export calculate_Phi_oe
export Double_hemisphere_Phi_oe
export Sph_with_Core_Phi_oe
export calculate_BSphere
export calculate_ASph_with_Core
export calculate_ESph_with_Core
export calculate_B
export calculate_C
export calculate_D
end # module SphrScatLib



