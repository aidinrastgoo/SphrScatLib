using LinearAlgebra
# in this test will set the value of SphrScatLib as normal Sphere with a PEC-Core and check if the result is the same.
ε_r1 = rand() # the both hiemsperes have same material 
ε_r2 = ε_r1 # + rand()
l_max = 100
radius_a = rand()
radius_b = radius_a + rand() #this radius shall be bigger as radius_a
Amplitude =rand() 
A2 = calculate_LH_A2(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max) # the coefficient for inside the SphrScatLib
A1 = calculate_HH_A1(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)
A_SC = calculate_ASph_with_Core(Amplitude,radius_a, radius_b, ε_r1, l_max) # the coefficient for inside Sphere with PEC_Core
#η = calculate_η(l_max, ε_r1, ε_r2)
#E_PEC = compute_E_PEC(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)
#E = calculate_OH_E(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max) # the coefficient for outside the SphrScatLib
#B = calculate_BSphere(Amplitude, radius_b , ε_r1, l_max, l_max) # the coofficient for the outside the normal Sphere
@test isapprox(A2' ,A_SC)   # if the coefficients are relaiable
@test isapprox(A1' ,A_SC)
#@test isapprox(E' , B)