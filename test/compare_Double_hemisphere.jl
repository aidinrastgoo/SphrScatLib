# in this test will set the value of SphrScatLib as Double hemisphere without a PEC-Core and check if the result is the same.
ε_r1 = 5 # the both hiemsperes have same material 
ε_r2 = 3
l_max = 4
radius_a = 0.0 #            ****there is no Metal-Sphere in the SphrScatLib****
radius_b = radius_a + 1.0    #this radius shall be bigger as radius_a
Amplitude =1.0 # rand() 
A2 = calculate_LH_A2(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max) # the coefficient for inside the SphrScatLib
A1 = calculate_HH_A1(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)
E = calculate_OH_E(Amplitude, radius_b, ε_r1, ε_r2, l_max)

B = calculate_B(Amplitude, radius_b, ε_r1, ε_r2,l_max)   #the coefficients of the potential outside the sphere ; radius = a
C = calculate_D(Amplitude, radius_b, ε_r1, ε_r2, l_max)   # the potential inside the double hemisphere
D = calculate_C(Amplitude, radius_b, ε_r1, ε_r2, l_max)  

@test isapprox(A2, C)    # if the coefficients are relaiable
@test isapprox(A1, D)
@test isapprox(E, B)