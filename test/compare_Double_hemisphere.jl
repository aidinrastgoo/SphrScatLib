# in this test will set the value of SphrScatLib as Double hemisphere without a PEC-Core and check if the result is the same.
using Plots
ε_r1 = rand(1:100) 
ε_r2 = rand(1:100)
l_max = 100
radius_a = 0.0 #            ****there is no Metal-Sphere in the SphrScatLib****
radius_b = radius_a + rand()    #this radius shall be bigger as radius_a
radius_out = radius_b + rand()
Amplitude = 1.0  
θ_degrees =  rand()
θ = θ_degrees * π / 180
ξ = cos(θ) 

A1 = calculate_HH_A1(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)[:,1] # the coefficient for region 1 of the SphrScatLib
C = calculate_C(Amplitude, radius_b, ε_r1, ε_r2, l_max)  # the coefficient for the northern hemisphere
A2 = calculate_LH_A2(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)[:,1] # the coefficient for region 2 of the SphrScatLib
D = calculate_D(Amplitude, radius_b, ε_r1, ε_r2, l_max)   # the coefficient for the southern hemisphere
@testset " Coefficient Comparison with Double hemisphere without PEC-Core" begin
    @test isapprox(A1, C)
    @test isapprox(A2, D)
end


r = zeros(Float64,l_max)
Phi = zeros(Float64,l_max)
Φ = zeros(Float64,l_max)
for i in 1:l_max-1
    r[i+1] = (i*0.05)
    Φ[i+1] =  Double_hemisphere_Phi(Amplitude, r[i+1] , radius_b, ξ, ε_r1, ε_r2, l_max)
    Phi[i+1] = calculate_Phi(Amplitude, r[i+1] , radius_a , radius_b, ξ ,ε_r1 ,ε_r1 , l_max)
end
labelfontsize = 5
labeltexts = ["SphrScatLib" "Double_hemisphere" ]
p = plot(r, [Phi Φ], 
         xlabel="Radius", ylabel="ϕ_1", label=labeltexts,
         title="Potential Comparison with Double_hemisphere",          titlefont=font(labelfontsize),          legendfont=font(labelfontsize),           guidefont=font(labelfontsize),            tickfont=font(labelfontsize),             textfont=font(labelfontsize), linewidth=[4 2])    


         vline!([r[1]], ylims=(minimum(Phi), Phi[1]), line=:dash, label="radius_a : $radius_a")
         vline!([radius_b], ylims=(minimum(Phi), Phi[1]), line=:dash, label="radius_b : $radius_b")
         vline!([r[10]], ylims=(minimum(Phi), Phi[1]), line=:dot, label="Amplitude : $Amplitude", linewidth=0.5, color=:white)
         vline!([r[1]], ylims=(minimum(Phi), Phi[1]), line=:dot, label="ε_r1 : $ε_r1 ", linewidth=0.5, color=:white)
         vline!([r[1]], ylims=(minimum(Phi), Phi[1]), line=:dot, label="ε_r2 : $ε_r2 ", linewidth=0.5, color=:white)
         vline!([r[1]], ylims=(minimum(Phi), Phi[1]), line=:dot, label="θ : $θ_degrees ", linewidth=0.5, color=:white)

 