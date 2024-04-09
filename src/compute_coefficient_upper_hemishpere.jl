

function calculate_HH_A1(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)
    A2 = calculate_LH_A2(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)
    A1 = zeros(Float64,l_max, 4)
    η= calculate_η(l_max, ε_r1, ε_r2)
    for i in 1: l_max
        A1[i, 1] = A2[i, 1] * η[i]
        A1[i, 2] = A2[i, 2]
        A1[i, 3] = A2[i, 3] * η[i]
        A1[1, 4] = A2[1, 4]
    end
    return A1
end