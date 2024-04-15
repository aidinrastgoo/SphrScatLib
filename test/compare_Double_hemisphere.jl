# in this test will set the value of SphrScatLib as Double hemisphere without a PEC-Core and check if the result is the same.
ε_r1 = rand() 
ε_r2 = rand()
l_max = 100
radius_a = 0.0 #            ****there is no Metal-Sphere in the SphrScatLib****
radius_b = radius_a + 1.0    #this radius shall be bigger as radius_a
radius_out = radius_b + rand()
Amplitude = 1.0 #rand() 
ξ = rand()

A1 = calculate_HH_A1(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)[:,1] # the coefficient for region 1 of the SphrScatLib
C = calculate_C(Amplitude, radius_b, ε_r1, ε_r2, l_max)  # the coefficient for the northern hemisphere
A2 = calculate_LH_A2(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)[:,1] # the coefficient for region 2 of the SphrScatLib
D = calculate_D(Amplitude, radius_b, ε_r1, ε_r2, l_max)   # the coefficient for the southern hemisphere


Φ_B = Double_hemisphere_Phi_oe(Amplitude, radius_out , radius_b, ξ, ε_r1, ε_r2, l_max)    # the potential outside the double hemisphere
Φ_E = calculate_Phi_oe(Amplitude, radius_out, radius_a , radius_b, ξ, ε_r1, ε_r2, l_max)  # the potential outside of the SphrScatLib




@testset " Coefficient Comparison with Double hemisphere without PEC-Core" begin
    @test isapprox(A1, C)
    @test isapprox(A2, D)
end

@testset " Potential of Enviroment Comparison with Double hemisphere without PEC-Core" begin
    @test Φ_B  ≈ Φ_E
end

