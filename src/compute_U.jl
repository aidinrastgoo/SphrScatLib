using SpecialFunctions
function A_kn(k, n)
    numerator_k = gamma(k/2 + 1) / gamma(k/2 + 1/2)
    denominator_n = gamma(n/2 + 1/2) / gamma(n/2 + 1)
    return numerator_k * denominator_n
end

function create_matrix(n, k)    # it's refixed matrix
    if n == k
        return 1.0 / (2n + 1)
    elseif (n + k) % 2 == 0
        return 0.0
    else
        t21 = (sin((π * n) / (2)) * cos((k * π) / (2))) / A_kn(k, n)
        t22 = n * (n + 1) - k * (k + 1)
        t23 = sin((k * π) / (2)) * cos((n * π) / (2)) * A_kn(k, n)
        return (2 / π) * ((t21 / t22) - (t23 / t22))
    end
end

function generate_U_matrix(l_max)
    U = zeros(Float64, l_max, l_max)
    for n in 0:l_max-1
        for k in 0:l_max-1
            U[n+1, k+1] = create_matrix(n, k)
        end
    end
    return U
end