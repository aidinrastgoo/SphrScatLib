using SphericalScattering
using StaticArrays
using LinearAlgebra
using SpecialFunctions  # For Legendre polynomials
using LegendrePolynomials
const μ0 = 4pi * 1e-7        # default permeability
const ε0 = 8.8541878176e-12  # default permittivity
# Define the direction of the incident electric field along the z-axis
dir = SVector(0.0, 0.0, 1.0) # Direction of the incident field along the z-axis ẑ = SVector(0.0, 0.0, 1.0)
# Define the amplitude of the incident electric field
amplitude = norm(dir) # Amplitude of the incident electric field
# Create the incident field (uniform electric field)
ex = UniformField(; direction=dir , amplitude)
# Define the parameters
r_a = 0.5 #(r_in + r_PEC) / 2   # Radius of the PEC core
r_b = 1.0 #(r_Out + r_in) / 2  # Radius of the outer sphere

# Define the surrounding medium (outer PEC-Core )
ε_r1 = ε_r2 = 5 #rand()    # Relative permittivity of the surrounding medium (assuming free space) ε_r1 = ε_surrounding / ε0
ε_surrounding = ε_r1 * ε0  
medium_surrounding = Medium(ε_surrounding, μ0)

sphere = LayeredSpherePEC(radii=SVector(r_b,r_a), filling = SVector(medium_surrounding))
# define observation points in both layers and outside of the sphere
r_in = (r_b + r_a) / 2 
r_out =  r_b * 2

θ = π / 4 # rand(0 : π)

OB_Die = SVector(0.0, r_in * sin(θ), r_in * cos(θ) )
OB_Out = SVector(0.0, r_out * sin(θ), r_out * cos(θ))

point_cart = [OB_Die , OB_Out] 



# Create the layered sphere with the PEC core
#sp = LayeredSphere(; radii=SVector(r_a, r_b), filling=SVector(medium_pec_core, medium_surrounding))
#sp = PECSphere(r_a)
sp = DielectricSphere(; radius= r_b, filling= medium_surrounding)

# the coefficient for inside and outside of Sphere with PEC_Core
l_max = 100

cos(θ)
# Compute the scalar potential at the observation point
Φ = scatteredfield(sphere, ex, ScalarPotential(point_cart))
Φ1 = scatteredfield(sp, ex, ScalarPotential(point_cart))
potential_in = calculate_Phi(amplitude, r_in , r_a, r_b, cos(θ) , ε_r1, ε_r2, l_max) 
potential_out = calculate_Phi(amplitude, r_out , r_a, r_b, cos(θ) , ε_r1, ε_r2, l_max) 
#@test Φ[1]  ≈ potential_in    # 
#@test (-Φ[1]  + Φ1[1])/2 ≈ potential_in
#@test Φ[2]  - ((amplitude * r_out) * cos(θ)) ≈ potential_out 



