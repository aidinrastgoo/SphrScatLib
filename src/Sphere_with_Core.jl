import .SphrScatLib: generate_σ
function calculate_MSph_Core(radius_a, radius_b, ε_r1, l_max)
    σ = generate_σ(l_max, l_max)
    M_S = zeros(Float64, l_max, l_max)
    
    for l in 0:l_max-1
        for m in 0:l_max-1
            term1 = (m+1) * ((radius_b^l) - ((radius_a^((2*l)+1))/(radius_b^(l+1))))
            term2 = ε_r1 * ((l* (radius_b^l)) + ((l +1)*((radius_a^((2*l)+1))/(radius_b^(l+1)))))
            M_S[l+1, m+1] = (2 / ((2*l) + 1)) * (term1 + term2) * σ[l+1, m+1]
        end
    end
    
    return M_S
end

function calculate_ASph(Amplitude, radius_b, l_max)
    σ = generate_σ(l_max, l_max)
    A_S = zeros(Float64, l_max)  # Adjusted dimensions to be a row vector
    for m in 0:l_max-1
        A_S[m+1] = - Amplitude * (((2 * (m + 2))/ 3 ) ) * radius_b * σ[2, m+1]
    end
    return A_S 
end
function calculate_ASph_with_Core(Amplitude, radius_a, radius_b, ε_r1, l_max)
    σ = generate_σ(l_max, l_max)
    A_S = calculate_ASph(Amplitude, radius_b, l_max)
    M_S = calculate_MSph_Core(radius_a, radius_b, ε_r1, l_max)
    return A_S / M_S
end
