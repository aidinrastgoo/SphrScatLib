using SphericalScattering
using Plots
using StaticArrays
using LinearAlgebra
const μ0 = 4pi * 1e-7        # default permeability
const ε0 = 8.8541878176e-12  # default permittivity
# Define the direction of the incident electric field along the z-axis
dir = SVector(0.0, 0.0, 1.0) # Direction of the incident field along the z-axis ẑ = SVector(0.0, 0.0, 1.0)

# Define the amplitude of the incident electric field
amplitude = norm(dir) # Amplitude of the incident electric field
# Create the incident field (uniform electric field)
ex = UniformField(; direction=dir , amplitude)
# Define the parameters
r_b = 1.0 # Radius of the outer sphere take as refrence radius and has value of 1.0 p.u 
r_a = round(rand(0.1:0.1:0.5), digits=1) #  # Radius of the PEC core is half of radius_a
radius = rand(0.5:0.1:2)#(r_a + r_b )/2
θ_degrees = rand(0:90)
θ = θ_degrees * π / 180
l_max = 40 

# Define the surrounding medium (outer PEC-Core )
ε_r1 = rand(1:99)             # Relative permittivity of the surrounding medium (assuming free space) ε_r1 = ε_surrounding / ε0
ε = ε_r1 * ε0  
Medium_1 = Medium(ε, μ0)
sphere = LayeredSpherePEC(radii=SVector(r_a,r_b), filling = SVector(Medium_1))

Φ = zeros(Float64,l_max)
r = zeros(Float64,l_max)
Phi = zeros(Float64,l_max)
x=0
for i in 1:l_max-1
    r[i+1] = (i*0.05)
    OB = SVector(0.0, r[i+1] * sin(θ), r[i+1] * cos(θ) ) # define observation points
    point_cart = [OB]  
    Φ1 = scatteredfield(sphere, ex, ScalarPotential(point_cart))
    Φ[i+1] = Φ1[1] - ((amplitude *r[i+1]) * cos(θ))
    Phi[i+1] = calculate_Phi(amplitude, r[i+1] , r_a, r_b, cos(θ) ,ε_r1 ,ε_r1 , l_max)
end
labelfontsize = 5
labeltexts = ["SphrScatLib" "SphericalScattering" ]
p = plot(r, [Phi Φ], 
         xlabel="Radius", ylabel="ϕ_1", label=labeltexts,
         title="Potential Comparison with SphericalScattering",          titlefont=font(labelfontsize),          legendfont=font(labelfontsize),           guidefont=font(labelfontsize),            tickfont=font(labelfontsize),             textfont=font(labelfontsize), linewidth=[4 2])    



         vline!([r_b], ylims=(minimum(Phi), Phi[31]), line=:dash, label="radius_b : $r_b")
         vline!([r_a], ylims=(minimum(Phi), Phi[10]), line=:dash, label="radius_a : $r_a")
         hline!([1], xlims=(r[1], r[end]), line=:dot, label="Amplitude : $amplitude", linewidth=0.5, color=:white)
         hline!([-amplitude], xlims=(r[1], r[end]), line=:dot, label="ε_r1 : $ε_r1 ", linewidth=0.5, color=:white)
         hline!([Phi[end]], xlims=(r[1], r[end]), line=:dot, label="θ : $θ_degrees ", linewidth=0.5, color=:white)

         
         
         
         

         
         


# Calculate Potential of southern hiemspere for each value of l
θ_degrees2 = rand(90:180)
θ2 = θ_degrees2 * π / 180
ε_r2 = rand(1:99)  
ϕ = zeros(l_max)
for l in 2:l_max+1
    ϕ[l-1] =  calculate_Phi(amplitude, radius , r_a, r_b, cos(θ2) ,ε_r1 ,ε_r2 , l)  # Stor elements
end
ϕ_value = round(ϕ[end], digits=6)

c = plot(1:l_max, [ϕ], xlabel="values of l", ylabel="ϕ_2", label="" , title="Potential Comparison in southern hiemspere for adjacent values of l",          titlefont=font(labelfontsize), legendfont=font(labelfontsize), guidefont=font(labelfontsize),   tickfont=font(labelfontsize),    textfont=font(labelfontsize)) 
hline!([ϕ[end]], xlims=(1, l_max), label="ϕ : $ϕ_value") 
hline!([minimum(ϕ)], xlims=(1, l_max), label="θ : $θ_degrees2", linewidth=0.5, color=:white) 
hline!([minimum(ϕ)], xlims=(1, l_max), label="r : $radius", linewidth=0.5, color=:white) 
hline!([minimum(ϕ)], xlims=(1, l_max), label="ε_r1 : $ε_r1", linewidth=0.5, color=:white) 
hline!([minimum(ϕ)], xlims=(1, l_max), label="ε_r2 : $ε_r2", linewidth=0.5, color=:white) 
hline!([minimum(ϕ)], xlims=(1, l_max), label="radius_a : $r_a", linewidth=0.5, color=:white) 
hline!([minimum(ϕ)], xlims=(1, l_max), label="radius_b : $r_b", linewidth=0.5, color=:white) 


display(plot(p, c))

@testset " Potential Comparison with SphericalScattering.jl" begin
    @test isapprox(Phi, Φ)
    println("Tolerance with SphericalScattering ", maximum(abs.(Phi - Φ)))
    if !isapprox(Phi, Φ)
        println("Tolerance exceeded for SphericalScattering by: ", maximum(abs.(A1 - D)))
    end
end