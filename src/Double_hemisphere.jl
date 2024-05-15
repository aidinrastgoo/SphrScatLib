
function calculate_B(Amplitude, radius, ε_r1, ε_r2,l_max)   #the coefficients of the potential outside the sphere ; radius = a
    U = generate_U_matrix(l_max)
    A = zeros(Float64, l_max)
    M = zeros(Float64, l_max, l_max)
    η = calculate_η(l_max, ε_r1, ε_r2)
    for k in 0:l_max-1
        for n in 0:l_max-1
           M[k+1, n+1] = (radius^(-(n + 2))) * ((η[k+1] * (n+1))  + (η[k+1] * k * ε_r1) + (((-1)^(n+k)) * (n+1)) + (((-1)^(n+k)) * k * ε_r2) ) * U[n+1, k+1]
           A[k+1] = Amplitude * ((η[k+1] * k * ε_r1) - η[k+1] + (((-1)^(1+k)) * k * ε_r2) - ((-1)^(1+k))) * U[2, k+1]
        end
    end
    return inv(M) * A
end


function calculate_D(Amplitude, radius, ε_r1, ε_r2, l_max)   # the potential inside the double hemisphere
    U = generate_U_matrix(l_max)
    A = zeros(Float64, l_max)
    M = zeros(Float64, l_max, l_max)
    η = calculate_η(l_max, ε_r1, ε_r2)
    for n in 0:l_max-1  
        for k in 0:l_max-1
            A[n+1] = -Amplitude * radius * ( (1/(n+1)) + 1 + (((-1)^(1+n)) / (n+1)) + ((-1)^(1+n)) ) * U[2, n+1]
            M[k+1, n+1] = radius^(k) * (η[k+1] + ((η[k+1] * k * ε_r1 )/(n+1)) +  ((-1)^(n+k)) +  ( (((-1)^(n+k)) * k * ε_r2)  /  (n+1))  ) * U[k+1, n+1] 
        end
    end
    return inv(M) * (A)  
end

function calculate_C(Amplitude, radius, ε_r1, ε_r2, l_max)
    D = calculate_D(Amplitude, radius, ε_r1, ε_r2, l_max)  
    η = calculate_η(l_max, ε_r1, ε_r2)
    return η.* D
end




function Double_hemisphere_Phi_oe(Amplitude, r , radius_b, ξ, ε_r1, ε_r2, l_max)
    Phi_oe = zeros(Float64, l_max)
    E = calculate_B(Amplitude, radius_b, ε_r2, ε_r1,l_max)
    for l in 0:l_max-1
        Phi_oe[l+1] =  ((E[l+1] * r^(-l-1))  * Pl(ξ,l) )  
    end
    return sum(Phi_oe) - (Amplitude* r * ξ)
end


    function Double_hemisphere_Phi_in(Amplitude, r , radius_b, ξ, ε_r1, ε_r2, l_max)
        if 0 <= ξ <= 1
            Phi_c = zeros(Float64, l_max)
            C = calculate_C(Amplitude, radius_b, ε_r1, ε_r2, l_max)
            #C = calculate_HH_A1(Amplitude, 0, radius_b, ε_r1, ε_r2, l_max)[:,1]
            for l in 0:l_max-1
                Phi_c[l+1] =  ((C[l+1] * r^l) )  * Pl(ξ,l) 
            end
            return  sum(Phi_c)
    
        elseif -1 <= ξ < 0
            Phi_d = zeros(Float64, l_max)
            D = calculate_D(Amplitude, radius_b, ε_r1, ε_r2, l_max)
           # D = calculate_LH_A2(Amplitude, 0, radius_b, ε_r1, ε_r2, l_max)[:,1]
            for l in 0:l_max-1
                Phi_d[l+1] =  ((D[l+1] * r^l))  * Pl(ξ,l) 
            end
            return  sum(Phi_d)
        else
            error("ξ has incorrect values outside the range")
        end
    end


function Double_hemisphere_Phi(Amplitude, r , radius_b, ξ, ε_r1, ε_r2, l_max)
    if r < radius_b
        return  Double_hemisphere_Phi_in(Amplitude, r , radius_b, ξ, ε_r1, ε_r2, l_max)  
    elseif r >= radius_b
        return  Double_hemisphere_Phi_oe(Amplitude, r , radius_b, ξ, ε_r1, ε_r2, l_max)
    else
        error("Review values: (Amplitude, r , rb, ξ, ε_r1, ε_r2, l_max)")
    end
end