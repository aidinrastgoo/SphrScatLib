include("compute_coefficient_lower_hemisphere.jl")
include("compute_η.jl")

function compute_E_PEC(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)
    A2 = calculate_LH_A2(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)
   # η = calculate_η(l_max, ε_r1, ε_r2)
    E = zeros(Float64, l_max)
    for l in 0:2:l_max-1

        if ε_r1 != ε_r2 
             E[l+1] = A2[l+1] * (radius_a^((2*l) + 1) ) * (l * ((l+1)^(-1))) 
                else E[l+1] = 0
        end
    
        
    end
    return E
end
