#inside
function calculate_ASph_with_Core(Amplitude, radius_a, radius_b, ε_r1, l_max)
    σ = I(l_max)
    A_S = zeros(Float64, l_max)
    M_S = zeros(Float64, l_max, l_max)
    for l in 0:l_max-1
        for k in 0:l_max-1
            A_S[k+1] = - Amplitude * radius_b * ((((k + 2) * 2)/ (((2*k) + 1) * (k+1))) )  * σ[2, k+1]
            term1 = ((radius_b^k) - ((radius_a^((2*k)+1))/(radius_b^(k+1))))
            term2 = (ε_r1 / (k+1)) * ((k* (radius_b^k)) + ((k + 1)*((radius_a^((2*k)+1))/(radius_b^(k+1)))))
            M_S[l+1, k+1] = (2 / ((2*l) + 1)) * (term1 + term2)  * σ[l+1, k+1]
        end
    end
#    A_S ./ M_S
    return inv(M_S) * A_S
end


# Outside


function calculate_ESph_with_Core(Amplitude, radius_b, ε_r1, l_max)
    σ = I(l_max)
    E_S = zeros(Float64, l_max)
    M_S_out = zeros(Float64, l_max, l_max)
    
    for l in 0:l_max-1
        for k in 0:l_max-1
            E_S[k+1] =  Amplitude * (ε_r1 - 1) * ( 2 /((2 * k) + 1)) * σ[2, k+1]
            M_S_out[l+1, k+1] = (2 / ((2*l) + 1)) * ((k * ε_r1) + k + 1) * (radius_b^(-k-2)) * σ[k+1, l+1]
        end
    end
return inv(M_S_out) * E_S
    #return E_S / M_S_out
end
