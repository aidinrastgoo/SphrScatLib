using LinearAlgebra
# in this test will set the value of SphrScatLib as normal Sphere and check if the result is the same.
ε_r2 = ε_r1 = .1 # the both hiemsperes have same material 
l_max = 100 
radius_a = 0.0   #            ****there is no Metal-Sphere in the SphrScatLib****
radius_b = radius_a + 1.0 #this radius shall be bigger as radius_a
Amplitude = rand() 
A2 = calculate_LH_A2(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max) # the coefficient for inside the SphrScatLib
A1 = calculate_HH_A1(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)
A_S = calculate_in(Amplitude, radius_b , ε_r1, l_max, l_max) # the coefficient for inside the normal Sphere
η = calculate_η(l_max, ε_r1, ε_r2)
E_PEC = compute_E_PEC(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)
E = calculate_OH_E(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max) # the coefficient for outside the SphrScatLib
B = calculate_BSphere(Amplitude, radius_b , ε_r1, l_max, l_max) # the coofficient for the outside the normal Sphere
@test isapprox(A2' ,A_S)   # if the coefficients are relaiable
@test isapprox(A1' ,A_S)
@test isapprox(B' ,E)