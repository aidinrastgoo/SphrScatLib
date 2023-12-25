using LinearAlgebra
ε_r2 = 1
ε_r1 = 1
l_max = 100
radius_a = 0.0
radius_b = radius_a + 0.9
Amplitude = 1.0
A2 = calculate_LH_A2(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)
A_S = calculate_in(Amplitude, radius, ε_r1, n_max, k_max)
@test isapprox(A2,A_S)



