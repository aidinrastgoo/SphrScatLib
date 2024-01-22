using LinearAlgebra
include("compute_U.jl")
include("compute_η.jl")
include("compute_E_PEC.jl")

function Part(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)
    η = calculate_η(l_max, ε_r1, ε_r2)
    U = generate_U_matrix(l_max)
    part = zeros(Float64, l_max, l_max)
    
    for l in 0:l_max-1  
        for m in 0:l_max-1   
            part[l+1, m+1] = -(((η[l+1])^2 * ε_r1 * (2*l + 1)) + (ε_r2 * ((-1)^(l+m)) * (2*l + 1))) * (radius_b^(-l-2)) * U[l+1, m+1] 
        end
    end 

    return part
end

function Part_PEC(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)
    E_PEC = compute_E_PEC(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)
    part = Part(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)

    return part * E_PEC 
end

function Part_Amplitude(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)
    η = calculate_η(l_max, ε_r1, ε_r2)
    U = generate_U_matrix(l_max)
    part = zeros(Float64, l_max)
    
    for l in 0:l_max-1  
        part[l+1] = Amplitude * ((η[l+1] * l * ε_r1) - η[l+1] + (ε_r2 * l * ((-1)^(l+1))) - ((-1)^(l+1))) * U[2, l+1] 
    end 

    return part
end

function Part_Out(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)
    η = calculate_η(l_max, ε_r1, ε_r2)
    U = generate_U_matrix(l_max)
    part = zeros(Float64, l_max, l_max)
    
    for m in 0:l_max-1  
        for l in 0:l_max-1   
            t1 = η[l+1] * ε_r1 * l 
            t2 = η[l+1] * (m+1)
            t3 = ε_r2 * l * ((-1)^(l+m))
            t4 = (m+1) * (-1)^(l+m)
            part[l+1, m+1] = (t1 + t2 + t3 + t4) * radius_b^(-m-2) * U[l+1, m+1]
        end
    end 

    return part
end

function calculate_OH_E(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)  #out of Hemisphere
    PEC = Part_PEC(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)
    Ampli = Part_Amplitude(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)
    Out = Part_Out(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)

    return inv(Out) * (Ampli + PEC)
end
