using LinearAlgebra
include("compute_U.jl")
include("compute_η.jl")
include("compute_E_PEC.jl")


function calculate_OH_E(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)  
    η = calculate_η(l_max, ε_r1, ε_r2)
    U = generate_U_matrix(l_max)
    term2 = zeros(Float64, l_max, l_max)
    term1 = zeros(Float64, l_max)
    for l in 0:l_max-1  
        for k in 0:l_max-1
            term1[k+1] = Amplitude  * ( (k * η[k+1] * ε_r1) - η[k+1] + (k * ε_r2 * ((-1)^(k+1))) - ((-1)^(k+1)) ) * U[2, k+1]
            term2[l+1 , k+1 ] =  ( (k * η[k+1] * ε_r1) + (η[k+1] * (l + 1)) + (k * ε_r2 * ((-1)^(l + k))) + ((l + 1) * ((-1)^(l + k))) ) * (radius_b^(-l-2)) * U[l+1, k+1] 
        end
    end 
term1
term2
    return inv(term2) * (term1)
end