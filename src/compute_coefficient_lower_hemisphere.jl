using LinearAlgebra
include("compute_U.jl")
include("compute_η.jl")

function calculate_LH_A2(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)  #lower Hiemspere
    η = calculate_η(l_max, ε_r1, ε_r2)
    U = generate_U_matrix(l_max)
    term2 = zeros(Float64, l_max, l_max)
    term1 = zeros(Float64, l_max)
    for m in 0:l_max-1  
        for l in 0:l_max-1
            term1[m+1] = - Amplitude * radius_b  * ( 1 + (1/(m+1)) + ( (-1)^(1+m) ) + (   (-1)^(1+m) * (1/(m+1))    )  ) * U[2, m+1]
            term2[l+1 , m+1 ] = ((radius_b^(l)) * ((η[l+1]  + ((l * η[l+1]  * ε_r1 ) / (m+1)) + ((-1)^(l+m)) + ((l * ε_r2 * ((-1)^(l+m))/(m+1))))) + ( (radius_b^(-l-1)) * (radius_a^((2 * l) + 1)) * (  ((l * η[l+1])/(l+1)) - ((l * η[l+1] * ε_r1)/(m+1)) + (l/(l+1)) * ((-1)^(l+m)) - ((l * ε_r2 * ((-1)^(l+m)))/(m+1))  ))) * U[l+1, m+1] 
        end
    end 
term1
term2
    return inv(term2) * (term1)
end