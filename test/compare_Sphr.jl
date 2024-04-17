
# in this test will set the value of SphrScatLib as normal Sphere and check if the result is the same.
ε_r2 = ε_r1 = 5 # the both hiemsperes have same material 
l_max = 100 
radius_a = 0.0   #            ****there is no Metal-Sphere****
radius_b = radius_a + 1.0 #this radius shall be bigger as radius_a
Amplitude = 1.0 
A2 = calculate_LH_A2(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)[:,1] # the coefficient for inside the SphrScatLib
A1 = calculate_HH_A1(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)[:,1]
A_S = calculate_in(Amplitude, radius_b , ε_r1, l_max, l_max)' # the coefficient for inside the normal Sphere

B = calculate_BSphere(Amplitude, radius_b , ε_r1, l_max, l_max)' # the coofficient for the outside of normal Sphere
E = calculate_HH_A1(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)[:,2] # the coofficient for the outside of SphrScatLib

@testset " Coefficient Comparison with Normal Sphere" begin
    @test isapprox(A1, A_S)
    @test isapprox(E, B)
    @test isapprox(A2, A_S)
end

