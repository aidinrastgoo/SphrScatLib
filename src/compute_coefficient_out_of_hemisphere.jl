using LinearAlgebra
include("compute_U.jl")
include("compute_η.jl")
include("compute_coefficient_lower_hemisphere.jl")

function calculate_first_Part(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)  
    η = calculate_η(l_max, ε_r1, ε_r2)
    U = generate_U_matrix(l_max)
    A2 = calculate_LH_A2(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)
    inside = zeros(Float64, l_max)
    effect = zeros(Float64, l_max)
    for l in 0:l_max-1  
        for m in 0:l_max-1
            # check the Sign if the result is not relaiable
            effect[l+1] =  Amplitude * U[l+1, 2] * ( radius_b * (  ( η[l+1] * ε_r1 * l) - η[l+1] + (ε_r2 * l * ((-1)^(l+1)))  - ((-1)^(l+1))  ))   
            inside[l+1] = U[l+1, m+1] * A2[l+1] * l * (radius_a^((2*l)+1)) * (radius_b^(-l-1)) * ( (η[l+1] * (1-ε_r1)) + (((-1)^(l+m)) * (1- ε_r2) ) )
        end
    end
    return effect - inside
end

function calculate_OH_E(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)  #out of Hiemspere
    η = calculate_η(l_max, ε_r1, ε_r2)
    U = generate_U_matrix(l_max)
    part1 = calculate_first_Part(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max) 
    part2 = zeros(Float64, l_max, l_max)
    
    for l in 0:l_max-1  
        for m in 0:l_max-1
            
            part2[l+1, m+1] = U[l+1, m+1] * (radius_b^(-l-1)) * ( (η[l+1] * ε_r1 * l)  + ( η[l+1] * (l+1) ) + ( ε_r2 * ((-1)^(l+m)) * l) + ( (l+1) * ((-1)^(l+m)) ) )
       
        end
    end
    return inv(part2) * part1
end
