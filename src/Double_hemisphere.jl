
function calculate_B(Amplitude, radius, ε_r1, ε_r2,l_max)   #the coefficients of the potential outside the sphere ; radius = a
    U = generate_U_matrix(l_max)
    A = zeros(Float64, l_max)
    M = zeros(Float64, l_max, l_max)
    η = calculate_η(l_max, ε_r1, ε_r2)
    for n in 0:l_max-1
        for k in 0:l_max-1
            term1 = η[k+1] * k * ε_r1
            term3 = ((-1)^(1 + k)) * k * ε_r2 
            term4 = (-1)^(1 + k)
            A[k+1] = Amplitude * (term1 - η[k+1] + term3 - term4) * U[2, k+1]
            t1 = η[k+1] * (n + 1)
            t2 = η[k+1] * k * ε_r1
            t3 = ((-1)^(n + k)) * (n + 1)
            t4 = ((-1)^(n + k)) * k * ε_r2
            M[n+1, k+1] = (radius^(-(n + 2))) * (t1 + t2 + t3 + t4) * U[n+1, k+1]
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