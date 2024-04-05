# in this test will set the value of SphrScatLib as Double hemisphere without a PEC-Core and check if the result is the same.
ε_r1 = rand() 
ε_r2 = rand()
l_max = 100
radius_a = 0.0 #            ****there is no Metal-Sphere in the SphrScatLib****
radius_b = radius_a + 1.0    #this radius shall be bigger as radius_a
Amplitude = rand() 
A2 = calculate_LH_A2(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)[:,1] # the coefficient for inside the SphrScatLib
A1 = calculate_HH_A1(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)[:,1] 
E = calculate_OH_E(Amplitude, radius_b, ε_r1, ε_r2, l_max)

B = calculate_B(Amplitude, radius_b, ε_r1, ε_r2,l_max)   #the coefficients of the potential outside the sphere ; radius = a
C = calculate_D(Amplitude, radius_b, ε_r1, ε_r2, l_max)   # the potential inside the double hemisphere
D = calculate_C(Amplitude, radius_b, ε_r1, ε_r2, l_max)  


@testset " Coefficient Comparison with Double hemisphere without PEC-Core" begin
    @test isapprox(A1, D)
    if !isapprox(A1, D)
        println("Tolerance exceeded for A1 and A_SC by: ", maximum(abs.(A1 - D)))
    end
    
    @test isapprox(E, B)
    if !isapprox(E, B)
        println("Tolerance exceeded for E and E_SC by: ", maximum(abs.(E - B)))
    end
    @test isapprox(A2, C)
    if !isapprox(A2, C)
        println("Tolerance exceeded by: ", maximum(abs.(A2 - C)))
    end

end

