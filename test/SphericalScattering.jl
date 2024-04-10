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
r_a = 0.5 #  # Radius of the PEC core is half of radius_a
radius = (r_a + r_b )/2
θ = π / 4 # rand(0 : π)
l_max = 40 

# Define the surrounding medium (outer PEC-Core )
ε_r1 = 3             # Relative permittivity of the surrounding medium (assuming free space) ε_r1 = ε_surrounding / ε0
ε = ε_r1 * ε0  
Medium_1 = Medium(ε, μ0)
sphere = LayeredSpherePEC(radii=SVector(r_b,r_a), filling = SVector(Medium_1))

Φ = zeros(Float64,l_max)
r = zeros(Float64,l_max)
Phi = zeros(Float64,l_max)

for i in 0:l_max-1
    r[i+1] = (i*0.05)
    OB = SVector(0.0, r[i+1] * sin(θ), r[i+1] * cos(θ) ) # define observation points
    point_cart = [OB]  
    Φ1 = scatteredfield(sphere, ex, ScalarPotential(point_cart))
    Φ[i+1] = Φ1[1] - ((amplitude *r[i+1]) * cos(θ))
    Phi[i+1] = calculate_Phi(amplitude, r[i+1] , r_a, r_b, cos(θ) ,ε_r1 ,ε_r1 , l_max)
end

labelfontsize = 5
labeltexts = ["SphrScatLib" "SphericalScattering"]
p = plot(r, [Phi Φ], 
         xlabel="Radius", ylabel="φ", label=labeltexts,
         title="Potential Comparison with SphericalScattering",          titlefont=font(labelfontsize),          legendfont=font(labelfontsize),           guidefont=font(labelfontsize),            tickfont=font(labelfontsize),             textfont=font(labelfontsize))    



# Calculate Potential of southern hiemspere for each value of l
ε_r2 = rand(1:99)  
ϕ = zeros(l_max)
for l in 2:l_max+1
    ϕ[l-1] =  calculate_Phi(amplitude, radius , r_a, r_b, cos(3π/4) ,ε_r1 ,ε_r2 , l)  # Stor elements
end

c = plot(1:l_max, [ϕ], xlabel="values of l", ylabel="ϕ", label="ε_r2 = $ε_r2", title="Potential Comparison in southern hiemspere for adjacent values of l",          titlefont=font(labelfontsize), legendfont=font(labelfontsize), guidefont=font(labelfontsize),   tickfont=font(labelfontsize),    textfont=font(labelfontsize))  



display(plot(p, c))