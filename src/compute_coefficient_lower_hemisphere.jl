

function calculate_LH_A2(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)  #lower Hiemspere
    η = calculate_η(l_max, ε_r1, ε_r2)
    U = generate_U_matrix(l_max)
    term2 = zeros(Float64, l_max, l_max)
    term1 = zeros(Float64, l_max)
    for l in 0:l_max-1  
        for k in 0:l_max-1
            term1[k+1] = - Amplitude * radius_b  * (  1 + (1/(k+1))  +  (-1)^(k+1) + (((-1)^(k+1)) / (k+1)) ) * U[2, k+1]
            term2[l+1 , k+1 ] =  (( ((η[l+1]) + (((-1)^(l+k)))) * ((radius_b^l) - ((radius_a^((2*l) + 1)) * radius_b^(-l-1) ))) + (  ((ε_r1 * η[l+1]) + (ε_r2 * (-1)^(l+k))) * ( ((l / (k+1))*radius_b^l) + (  ((l+1)/(k+1)) * (radius_a^((2*l)+1)) * (radius_b^(-l-1)) )))) * U[l+1, k+1] 
        end
    end 
    return inv(term2) * (term1)
end