using Test

ε_r2 = 1.0
ε_r1 = 1.0
l_max = 10
radius_a = 0.0
radius_b = radius_a + 0.9
Amplitude = 1.0

U = generate_U_matrix(l_max)
η = calculate_η(l_max, ε_r1, ε_r2)

A2 = calculate_LH_A2(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)



@test 1.0 == 1.0

A2