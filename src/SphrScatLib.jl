module SphrScatLib

import Pkg
include("compute_η.jl")
include("compute_U.jl")
include("compute_coefficient_lower_hemisphere.jl")
include("compute_coefficient_out_of_hemisphere.jl")
include("Sphere.jl")


export calculate_η
export generate_U_matrix
export calculate_LH_A2
export calculate_OH_E
export calculate_in

end # module SphrScatLib
