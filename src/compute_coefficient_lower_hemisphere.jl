
function coefficients(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)  #lower Hiemspere
    η = calculate_η(l_max, ε_r1, ε_r2)
    U = generate_U_matrix(l_max)
    if radius_a == 0.0                          # absence of a core
        B = zeros(Float64, 1 , 2*l_max)                   
        A =  zeros(Float64, 2*l_max , 2*l_max)    #eliminating the 2nd term due to the absence of a core
        C = zeros(Float64,   2 * l_max, 1)        #eliminating coefficients E_PEC & A_PEC
        for k in 0:l_max-1  
            for l in 0 :l_max-1
                B[k+1] = -Amplitude * radius_b * ( (1/(k+1)) + 1 + (((-1)^(1+k)) / (k+1)) + ((-1)^(1+k)) ) * U[2, k+1]  # first response vector
                A[k+1,       l+1 ] = ((radius_b^(l)) *  (   η[l+1] + ((-1)^(l+k)) +  ((l * ε_r1 * η[l+1] ) / (k+1)) + ( (l * ε_r2 * ((-1)^(l+k)))  / (k+1) ) ) )* U[k+1, l+1]                # 1st term   
                B[l_max+k+1] = Amplitude  * ( (k * η[k+1] * ε_r1) - η[k+1] + (k * ε_r2 * ((-1)^(k+1))) - ((-1)^(k+1)) ) * U[2, k+1] # second response vector
                A[l_max+k+1 , l_max+l+1 ] =   ( (k * η[k+1] * ε_r1) + (η[k+1] * (l + 1)) + (k * ε_r2 * ((-1)^(l + k))) + ((l + 1) * ((-1)^(l + k))) ) * (radius_b^(-l-2)) * U[k+1, l+1]  # 3rd term
            end
        end 
        C = (B / A)   #coefficient
        return   C

        elseif radius_a > 0.0
            B = zeros(Float64, (3*l_max)+1)
            A = zeros(Float64, (3*l_max)+1 , (3*l_max)+1)
            for k in 0:l_max-1  
                for l in 0 :l_max-1
                    B[k+1] = -Amplitude * radius_b * ( (1/(k+1)) + 1 + (((-1)^(1+k)) / (k+1)) + ((-1)^(1+k)) ) * U[ 2 ,k+1]  # first response vector
                    A[k+1,       l+1 ] =  ((radius_b^(l)) *  (   η[l+1] + ((-1)^(l+k)) +  ((l * ε_r1 * η[l+1] ) / (k+1)) + ( (l * ε_r2 * ((-1)^(l+k)))  / (k+1) ) ) )* U[k+1, l+1]                 # 1st term 
                    A[k+1,       l+(2*l_max)+1] =  (radius_b^(-l-1)) * (  η[l+1] + ((-1)^(l+k)) -  (((l+1) * ε_r1 * η[l+1] ) / (k+1)) - ( ((l+1) * ε_r2 * ((-1)^(l+k)))  / (k+1) ) )* U[k+1, l+1]     # 2nd term
                    
                    B[l_max+k+1] = Amplitude  * ( (k * η[k+1] * ε_r1) - η[k+1] + (k * ε_r2 * ((-1)^(k+1))) - ((-1)^(k+1)) ) * U[2 , k+1]   # second response vector
                    A[k+l_max+1 , l+(l_max)+1 ] =  ( (k * η[k+1] * ε_r1) + (η[k+1] * (l + 1)) + (k * ε_r2 * ((-1)^(l + k))) + ((l + 1) * ((-1)^(l + k))) ) * (radius_b^(-l-2)) * U[k+1, l+1]  # 3rd term
                    A[k+l_max+1, l+(2*l_max)+1] =   -(((l+1) * η[k+1] * η[l+1] * ε_r1) + (k * η[k+1] * η[l+1] * ε_r1) + (k * ε_r2 * ((-1)^(l + k))) + ((l + 1) * ε_r2  * ((-1)^(l + k)) )) * (radius_b^(-l-2)) * U[k+1, l+1] # 4th term
                                    
                    A[k+(2*l_max)+1, l+1]  = (radius_a^l) * ((η[l+1]) + ((-1)^(l+k)) ) * U[k+1, l+1]    # 5th term
                    A[k+(2*l_max)+1, l+(2*l_max)+1] = (radius_a^(-l-1)) * ((η[l+1]) + ((-1)^(l+k)) ) * U[k+1, l+1]  # 6th term
                    A[k+(2*l_max)+1, end] =  - (1 + (-1)^k) * U[k+1, 1]     # 7th term
                    A[(3*l_max)+1,       (2*l_max)+1] = 1    # nullifying E_PEC in l =  0 
                end
            end 
            C = (A \ B)
            return  C
        end
end

function calculate_LH_A2(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)  # Arranging coefficients of lower Hiemspere
    co = coefficients(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)
    x = zeros(Float64, l_max , 4)
    for i in 1:l_max
        x[i, 1] = co[i]
        x[i, 2] = co[i+l_max]
        if size(co, 1) == (3 * l_max) + 1
            x[i, 3] = co[i + (2*l_max)]
            x[1, 4] = co[end]
        end
    end
    return x
end


function coefficients_EDL(Amplitude, radius_a, radius_b, ε_r1, ε_r2, ε_m, Δ, l_max)  #impedance layer
    η = calculate_η(l_max, ε_r1, ε_r2)
    U = generate_U_matrix(l_max)
    if radius_a == 0.0
        error("PEC size has incorrect value")
    else
        B = zeros(Float64, (3*l_max)+1)
        A = zeros(Float64, (3*l_max)+1 , (3*l_max)+1)
        for k in 0:l_max-1  
            for l in 0 :l_max-1
                B[k+1] = -Amplitude * radius_b * ( (1/(k+1)) + 1 + (((-1)^(1+k)) / (k+1)) + ((-1)^(1+k)) ) * U[ 2 ,k+1]  # first response vector
                A[k+1,       l+1 ] =  ((radius_b^(l)) *  (   η[l+1] + ((-1)^(l+k)) +  ((l * ε_r1 * η[l+1] ) / (k+1)) + ( (l * ε_r2 * ((-1)^(l+k)))  / (k+1) ) ) )* U[k+1, l+1]                 # 1st term 
                A[k+1,       l+(2*l_max)+1] =  (radius_b^(-l-1)) * (  η[l+1] + ((-1)^(l+k)) -  (((l+1) * ε_r1 * η[l+1] ) / (k+1)) - ( ((l+1) * ε_r2 * ((-1)^(l+k)))  / (k+1) ) )* U[k+1, l+1]     # 2nd term
                
                B[l_max+k+1] = Amplitude  * ( (k * η[k+1] * ε_r1) - η[k+1] + (k * ε_r2 * ((-1)^(k+1))) - ((-1)^(k+1)) ) * U[2 , k+1]   # second response vector
                A[k+l_max+1 , l+(l_max)+1 ] =  ( (k * η[k+1] * ε_r1) + (η[k+1] * (l + 1)) + (k * ε_r2 * ((-1)^(l + k))) + ((l + 1) * ((-1)^(l + k))) ) * (radius_b^(-l-2)) * U[k+1, l+1]  # 3rd term
                A[k+l_max+1, l+(2*l_max)+1] =   -(((l+1) * η[k+1] * η[l+1] * ε_r1) + (k * η[k+1] * η[l+1] * ε_r1) + (k * ε_r2 * ((-1)^(l + k))) + ((l + 1) * ε_r2  * ((-1)^(l + k)) )) * (radius_b^(-l-2)) * U[k+1, l+1] # 4th term
                                
                A[k+(2*l_max)+1, l+1]  = ((Δ * ((((ε_r1 * η[l+1]) + (ε_r2 * ((-1)^(k+l)))) * (l) * radius_a^(l-1)) / (ε_m)))             -((radius_a^l) * ((η[l+1]) + ((-1)^(l+k)) )) ) * U[k+1, l+1]  # 5th term
                A[k+(2*l_max)+1, l+(2*l_max)+1] = - ((Δ * ((((ε_r1 * η[l+1]) + (ε_r2 * ((-1)^(k+l)))) * (l+1) * radius_a^(-l-2)) / (ε_m)))            + ( (radius_a^(-l-1)) * ((η[l+1]) + ((-1)^(l+k)) ))) * U[k+1, l+1]  # 6th term
                A[k+(2*l_max)+1, end] =  (1 + (-1)^k) * U[k+1, 1]     # 7th term
                A[(3*l_max)+1,       (2*l_max)+1] = 1    # nullifying E_PEC in l =  0 
            end
        end    
        C = A \ B
        return C
    end
end




function calculate_EDL_A2(Amplitude, radius_a, radius_b, ε_r1, ε_r2, ε_m, Δ, l_max)  # Arranging coefficients in an impedance layer example
    co = coefficients_EDL(Amplitude, radius_a, radius_b, ε_r1, ε_r2, ε_m, Δ, l_max)
    x = zeros(Float64, l_max , 4)
    for i in 1:l_max
        x[i, 1] = co[i]
        x[i, 2] = co[i+l_max]
        if size(co, 1) == (3 * l_max) + 1
            x[i, 3] = co[i + (2*l_max)]
            x[1, 4] = co[end]
        end
    end
    return x
end
