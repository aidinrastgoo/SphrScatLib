using LinearAlgebra
# in this test will set the value of SphrScatLib as normal Sphere and check if the result is the same.
ε_r2 = ε_r1 = rand() # the both hiemsperes have same material 
l_max = 100 
radius_a = 0.0    #            ****there is no Metal-Sphere in the SphrScatLib****
radius_b = radius_a + rand() #this radius shall be bigger as radius_a
Amplitude = rand() 
A2 = calculate_LH_A2(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max) # the coofficient for inside the SphrScatLib
A_S = calculate_in(Amplitude, radius_b , ε_r1, l_max, l_max) # the coofficient for inside the normal Sphere
E = calculate_OH_E(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max) # the coofficient for outside the SphrScatLib
B = calculate_in(Amplitude, radius_b , ε_r1, l_max, l_max) # the coofficient for the outside the normal Sphere
@test isapprox(A2' ,A_S) || isapprox(B' ,E) # if the coofficients are relaiable