using SpecialFunctions
function A_lk(l, k)
    numerator_l = gamma(l/2 + 1) / gamma(l/2 + 1/2)
    denominator_k = gamma(k/2 + 1/2) / gamma(k/2 + 1)
    return numerator_l * denominator_k
end

function create_matrix(l, k)    # it's refixed matrix
    if l == k
        return 1.0 / (2l + 1)
    elseif (l + k) % 2 == 0
        return 0.0
    else
        return ((A_lk(l, k) * sin((k * π) / (2)) * cos((l * π) / (2))) - ((1/ A_lk(l, k)) * sin((l * π) / (2)) * cos((k * π) / (2))) ) / ((π/2) * (k-l) * (k + l +1))
    end
end

function generate_U_matrix(l_max)
    U = zeros(Float64, l_max, l_max)
    for l in 0:l_max-1
        for k in 0:l_max-1
            U[l+1, k+1] = create_matrix(l, k)
        end
    end
    return U
end