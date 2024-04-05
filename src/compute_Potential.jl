

function calculate_Phi_1(Amplitude, r , ra, rb, ξ, ε_r1, ε_r2, l_max)
    Phi_1 = zeros(Float64, l_max)
    coefficients = calculate_HH_A1(Amplitude, ra, rb, ε_r1, ε_r2, l_max)
    Al = coefficients[:,1]
    E = coefficients[:,2]
    for l in 0:l_max-1
        Phi_1[l+1] =  ((Al[l+1] * r^l) + (E[l+1, 1] * r^(-(l+1))))  * Pl(ξ,l) 
    end
    return sum(Phi_1)
end

function calculate_Phi_2(Amplitude, r , ra, rb, ξ, ε_r1, ε_r2, l_max)
    Phi_1 = zeros(Float64, l_max)
    coefficients = calculate_LH_A2(Amplitude, ra, rb, ε_r1, ε_r2, l_max)
    Al = coefficients[:,1]
    E = coefficients[:,2]
    for l in 0:l_max-1
        Phi_1[l+1] =  ((Al[l+1] * r^l) + (E[l+1, 1] * r^(-(l+1))))  * Pl(ξ,l) 
    end
    return sum(Phi_1)
end

function calculate_Phi_o(Amplitude, r , rb, ξ, ε_r1, ε_r2, l_max)
    Phi_o = zeros(Float64, l_max)
    E = calculate_OH_E(Amplitude,  rb, ε_r1, ε_r2, l_max)
    for l in 0:l_max-1
        Phi_o[l+1] =  ((E[l+1] * r^(-l-1)))  * Pl(ξ,l) 
    end
    return sum(Phi_o)
end

function calculate_Phi_oe(Amplitude, r , rb, ξ, ε_r1, ε_r2, l_max)
    Phi_oe = zeros(Float64, l_max)
    E = calculate_OH_E(Amplitude,  rb, ε_r1, ε_r2, l_max)
    for l in 0:l_max-1
        Phi_oe[l+1] =  ((E[l+1] * r^(-l-1))  * Pl(ξ,l) ) - (Amplitude * r * Pl(ξ,1)) 
    end
    return sum(Phi_oe)
end

function calculate_Phi_PEC(Amplitude, ra ,  rb , ε_r1, ε_r2, l_max)
    coefficients = calculate_HH_A1(Amplitude, ra, rb, ε_r1, ε_r2, l_max)
    return coefficients[3,1]
end


function calculate_Phi(Amplitude, r, ra, rb, ξ, ε_r1, ε_r2, l_max)
    if r <= ra
        return calculate_Phi_PEC(Amplitude, ra, rb, ε_r1, ε_r2, l_max)
    elseif ra < r < rb && 0 <= ξ <= 1
        return calculate_Phi_1(Amplitude, r, ra, rb, ξ, ε_r1, ε_r2, l_max)
    elseif ra < r < rb && 0 < ξ <= -1
        return calculate_Phi_2(Amplitude, r, ra, rb, ξ, ε_r1, ε_r2, l_max)
    elseif r >= rb
        return calculate_Phi_oe(Amplitude, r, rb, ξ, ε_r1, ε_r2, l_max)
    else
        error("Review values: (Amplitude, r, ra, rb, ξ, ε_r1, ε_r2, l_max)")
    end
end



        



