function create_σ(n, k)
    if n == k
        return 1
    else
        return 0
    end
end

function generate_σ(n_max, k_max)
    σ = zeros(Int, n_max, k_max)
    for n in 1:n_max
        for k in 1:k_max
            σ[n, k] = create_σ(n, k)
        end
    end
    return σ
end

function calculate_MSph(radius, ε_r1, n_max, k_max)
    σ = generate_σ(n_max, k_max)
    M_s = zeros(Float64, n_max, k_max)
    
    for n in 0:n_max-1
        for k in 0:k_max-1
            term1 = - n * ε_r1 
            term2 = -(n + 1)
            M_s[n+1, k+1] = (radius^(-(n + 2))) * (term1 + term2) * σ[n+1, k+1]
        end
    end
    
    return M_s
end

function calculate_ASph(Amplitude, ε_r1, k_max)
    σ = generate_σ(k_max, k_max)
    A_s = zeros(Float64, 1, k_max)  # Adjusted dimensions to be a row vector
    for k in 0:k_max-1
        A_s[1, k+1] = Amplitude * (1-(k * ε_r1)) * σ[2, k+1]
    end
    
    return A_s
end

function calculate_BSphere(Amplitude, radius, ε_r1, n_max, k_max)
    M_s = calculate_MSph(radius, ε_r1, n_max, k_max)
    A_s = calculate_ASph(Amplitude, ε_r1, k_max)
    B_s = A_s/ M_s
    return B_s
end

function calculate_in(Amplitude, radius, ε_r1, n_max, k_max)
    B = calculate_BSphere(Amplitude, radius, ε_r1, n_max, k_max)
    σ = generate_σ(k_max, k_max)

    S1 = zeros(Float64, k_max)
    S2 = zeros(Float64, k_max)
    S3 = zeros(Float64, k_max, k_max)

    for k in 0:k_max-1
        S1[k+1] = B[1, k+1] * radius^(-k-1)  # Accessing B as a matrix
        S2[k+1] = radius * Amplitude * σ[2, k+1]
        S3[k+1, k+1] = (radius^k) * σ[k+1, k+1]
    end

    S = (S1' * σ - S2') / S3  # Matrix division using /
    return S
end