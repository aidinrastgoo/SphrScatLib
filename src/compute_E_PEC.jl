include("compute_coefficient_lower_hemisphere.jl")
include("compute_η.jl")
include("compute_U.jl")

function compute_E_PEC(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)
    A2 = calculate_LH_A2(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)
    U = generate_U_matrix(l_max)
    part_A2 = zeros(Float64, l_max)

    for l in 0:l_max-1
        
        part_A2[l+1] = -A2[l+1] * ((l * (radius_a^((2*l) +1)))/(l+1)) 
        
    end
    return  part_A2 
end
