
# in this test will set the value of SphrScatLib as normal Sphere with a PEC-Core and check if the result is the same.
ε_r1 = rand() 
ε_r2 = ε_r1  ## the both hiemsperes have same material 
ε_m = rand()
l_max = 50
Δ = 0 # If  approach to zero results in the same coefficient
radius_a = rand()
radius_b = radius_a + 1.0 #this radius shall be bigger as radius_a
radius_out = radius_b + rand()
ξ = rand()
Amplitude =rand() 
A2 = calculate_LH_A2(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)[:,1] # the coefficient for inside the SphrScatLib
A1 = calculate_HH_A1(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)[:,1]
A2_EDL = calculate_EDL_A2(Amplitude, radius_a, radius_b, ε_r1, ε_r2, ε_m, Δ, l_max)[:,1]
Φ_E = calculate_Phi_oe(Amplitude, radius_out, radius_a , radius_b, ξ, ε_r1, ε_r2, l_max) 


A_SC = calculate_ASph_with_Core(Amplitude,radius_a, radius_b, ε_r1, l_max) # the coefficient for inside Sphere with PEC_Core
Φ_SC = Sph_with_Core_Phi_oe(Amplitude, radius_out ,radius_a , radius_b, ξ, ε_r1)

@testset " Coefficient Comparison with Sphere with a PEC-Core" begin
    @test isapprox(A1, A_SC)
    @test isapprox(A2, A_SC)
    @test isapprox(A2_EDL, A_SC)
end

@testset " Enviroment Potential Comparison with Sphere with a PEC-Core" begin
    @test isapprox(Φ_SC, Φ_E)
end

