
# in this test will set the value of SphrScatLib as normal Sphere and check if the result is the same.
ε_r2 = ε_r1 = 5 # the both hiemsperes have same material 
l_max = 100 
radius_a = 0.0   #            ****there is no Metal-Sphere in the SphrScatLib****
radius_b = radius_a + 1.0 #this radius shall be bigger as radius_a
Amplitude = 1.0 
A2 = calculate_LH_A2(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)[:,1] # the coefficient for inside the SphrScatLib
A1 = calculate_HH_A1(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)[:,1]
A_S = calculate_in(Amplitude, radius_b , ε_r1, l_max, l_max)' # the coefficient for inside the normal Sphere
η = calculate_η(l_max, ε_r1, ε_r2)
E = calculate_OH_E(Amplitude, radius_b, ε_r1, ε_r2, l_max) # the coefficient for outside the SphrScatLib
B = calculate_BSphere(Amplitude, radius_b , ε_r1, l_max, l_max)' # the coofficient for the outside the normal Sphere


@testset " Coefficient Comparison with Normal Sphere" begin
    @test isapprox(A1, A_S)
    if !isapprox(A1, A_S)
        println("Tolerance exceeded for A1 and A_SC by: ", maximum(abs.(A1 - A_S)))
    end
    
    @test isapprox(E, B)
    if !isapprox(E, B)
        println("Tolerance exceeded for E and E_SC by: ", maximum(abs.(E - B)))
    end
    @test isapprox(A2, A_S)
    if !isapprox(A2, A_S)
        println("Tolerance exceeded by: ", maximum(abs.(A2 - A_S)))
    end

end

