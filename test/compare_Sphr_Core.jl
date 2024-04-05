
# in this test will set the value of SphrScatLib as normal Sphere with a PEC-Core and check if the result is the same.
ε_r1 = rand() 
ε_r2 = ε_r1  ## the both hiemsperes have same material 
l_max = 100
radius_a = rand()
radius_b = radius_a + 1.0 #this radius shall be bigger as radius_a
Amplitude =rand() 
A2 = calculate_LH_A2(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)[:,1] # the coefficient for inside the SphrScatLib
A1 = calculate_HH_A1(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)[:,1]
E = calculate_OH_E(Amplitude, radius_b, ε_r1, ε_r2, l_max)
A_SC = calculate_ASph_with_Core(Amplitude,radius_a, radius_b, ε_r1, l_max) # the coefficient for inside Sphere with PEC_Core
E_SC = calculate_ESph_with_Core(Amplitude, radius_b, ε_r1, l_max)

tolerance = 1e-2  # Adjust tolerance as needed

@testset " Coefficient Comparison with Sphere with Core" begin
    @test isapprox(A1, A_SC, atol=tolerance)
    if !isapprox(A1, A_SC, atol=tolerance)
        println("Tolerance exceeded for A1 and A_SC by: ", maximum(abs.(A1 - A_SC)))
    end
    
    @test isapprox(E, E_SC, atol=tolerance)
    if !isapprox(E, E_SC, atol=tolerance)
        println("Tolerance exceeded for E and E_SC by: ", maximum(abs.(E - E_SC)))
    end
    @test isapprox(A2, A_SC, atol=tolerance)
    if !isapprox(A2, A_SC, atol=tolerance)
        println("Tolerance exceeded by: ", maximum(abs.(A2 - A_SC)))
    end

end