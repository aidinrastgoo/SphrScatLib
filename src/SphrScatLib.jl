module SphrScatLib

import Pkg
include("compute_η.jl")
include("compute_U.jl")
include("compute_coefficient_lower_hemisphere.jl")


export calculate_η
export generate_U_matrix
export calculate_LH_A2

end # module SphrScatLib
