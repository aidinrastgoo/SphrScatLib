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


function Sph_with_Core_Phi_oe(Amplitude, r ,radius_a , radius_b, ξ, ε_r1)
    
    E = (Amplitude/((ε_r1+2) * r^2)) * (    ((9 * ε_r1) / ( ( (2 /radius_a^3) - (2/radius_b^3) ) + ( ( ε_r1 / radius_a^3 ) + ((2 * ε_r1 )/ (radius_b^3) ) )  )  )   + (  (ε_r1 - 1)  * (radius_b^3 )  )     - ((ε_r1+2)r^3   )       )  * ξ
    return E
end
