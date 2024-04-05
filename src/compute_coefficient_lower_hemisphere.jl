using IterativeSolvers

function coefficients(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)  #lower Hiemspere
    η = calculate_η(l_max, ε_r1, ε_r2)
    U = generate_U_matrix(l_max)
    if radius_a == 0.0
        term1 = zeros(Float64,1 ,  l_max)
        term2 = zeros(Float64, l_max , l_max)
        for k in 0:l_max-1  
            for l in 0 :l_max-1
                term1[1, k+1] = -Amplitude * radius_b * ( (1/(k+1)) + 1 + (((-1)^(1+k)) / (k+1)) + ((-1)^(1+k)) ) * U[2, k+1] 
                term2[k+1,       l+1 ] =  ((radius_b^(l)) *  (   η[l+1] + ((-1)^(l+k)) +  ((l * ε_r1 * η[l+1] ) / (k+1)) + ( (l * ε_r2 * ((-1)^(l+k)))  / (k+1) ) ) )* U[k+1, l+1]                 # A_l (2)     r_b
            end
        end 
        x = term1 / term2
        return x'
        elseif radius_a > 0.0
            term1 = zeros(Float64,1 ,  (2*l_max)+1)
            term2 = zeros(Float64, (2*l_max)+1 , (2*l_max)+1)
            for k in 0:l_max-1  
                for l in 0 :l_max-1
                    term1[1, k+1] = -Amplitude * radius_b * ( (1/(k+1)) + 1 + (((-1)^(1+k)) / (k+1)) + ((-1)^(1+k)) ) * U[2, k+1] 
                    term2[k+1,       l+1 ] =  ((radius_b^(l)) *  (   η[l+1] + ((-1)^(l+k)) +  ((l * ε_r1 * η[l+1] ) / (k+1)) + ( (l * ε_r2 * ((-1)^(l+k)))  / (k+1) ) ) )* U[k+1, l+1]                 # A_l (2)     r_b
                    term2[k+1,       l+l_max+1] =  (radius_b^(-l-1)) * (  η[l+1] + ((-1)^(l+k)) -  (((l+1) * ε_r1 * η[l+1] ) / (k+1)) - ( ((l+1) * ε_r2 * ((-1)^(l+k)))  / (k+1) ) )* U[k+1, l+1]     # E_l (PEC)    r_b
                    term2[k+l_max+1, l+1]  = (radius_a^l) * ((η[l+1]) + ((-1)^(l+k)) ) * U[k+1, l+1]                                                                           # A_l (2)     r_a
                    term2[k+l_max+1, l+l_max+1] = (radius_a^(-l-1)) * ((η[l+1]) + ((-1)^(l+k)) ) * U[k+1, l+1]                                                                 #E_l (PEC)    r_a
                    term2[k+l_max+1, end] =  - (1 + (-1)^k) * U[k+1, 1]
                    term2[end,       l_max+1] = 1
                end
            end 
            x = IterativeSolvers.gmres(term2,term1;verbose=true, log=true)[1]
            return  IterativeSolvers.gmres(term2,term1;verbose=true, log=true)[1]
        end
end

function calculate_LH_A2(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)
    co = coefficients(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)
    x = zeros(Float64, l_max , 3)
    for i in 1:l_max
        x[i, 1] = co[i]
        if size(co, 1) == (2 * l_max) + 1
            x[i, 2] = co[i + l_max]
            x[1, 3] = co[end]
        end
    end
    return x
end