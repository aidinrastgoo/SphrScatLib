using LinearAlgebra
include("compute_U.jl")
include("compute_η.jl")
function Parts(radius_a, radius_b, ε_r1, ε_r2, l_max)
    η = calculate_η( l_max, ε_r1, ε_r2)
    U = generate_U_matrix(l_max)
    part1 = zeros(Float64, l_max, l_max)
    for l in 0:l_max-1  
        for m in 0:l_max-1
            if ε_r1 == ε_r2 
                part1[l+1 , m+1 ] = U[l+1, m+1] * ((((η[l+1] * (l+1) ) + ( ((-1)^(l+m)) * (l + 1) ) + ( η[l+1] * ε_r1 * l ) + ( ((-1)^(l+m)) * ε_r2 * l)) * (radius_b^l))+(((radius_a^((2 * l) + 1)) * l * ( (1 - ε_r1) + ((1 - ε_r2) * (-1)^(l+m)) ) ) * (radius_b^(-l-1))))
            else part1[l+1 , m+1 ] = U[l+1, m+1] * ((((η[l+1] * (l+1) ) + ( ((-1)^(l+m)) * (l + 1) ) + ( η[l+1] * ε_r1 * l ) + ( ((-1)^(l+m)) * ε_r2 * l)) * (radius_b^l))+(((radius_a^((2 * l) + 1)) * l * ( (1 - ε_r1) + ((1 - ε_r2) * (-1)^(l+m)) ) ) * (radius_b^(-l-1)) * (((η[l+1] * ε_r1 ) - (ε_r2 * (-1)^(l+m))) * ((ε_r1 ) - (ε_r2 * (-1)^(l+m)))^(-1))))
            end
        end
    end 
    return part1

end

function calculate_LH_A2(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)  #lower Hiemspere
    η = calculate_η(l_max, ε_r1, ε_r2)
    U = generate_U_matrix(l_max)
    term1 = zeros(Float64, l_max)
    term2 = Parts(radius_a, radius_b, ε_r1, ε_r2, l_max)
    for l in 0:l_max-1  
        for m in 0:l_max-1
            term1[l+1] = - Amplitude * U[l+1, 2] * ( radius_b * (1 + ( (-1)^(1+l) ) + (  ((-1)^(1+l)) * (l+1) )  + (l+1)     ))

           # term2[l+1, m+1] = U[l+1, m+1] * (((radius_b^l) * part1[l+1 , m+1]) +   ((radius_b^(-l-1)) * part2[l+1 + m+1 ]) )
        end
    end
    return inv(term2) * (term1)
end