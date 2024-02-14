using SphericalScattering
using StaticArrays
using LinearAlgebra
const μ0 = 4pi * 1e-7        # default permeability
const ε0 = 8.8541878176e-12  # default permittivity
# Define the direction of the incident electric field along the z-axis
dir = SVector(1.0, 0.5, -3.0) # Direction of the incident field along the z-axis ẑ = SVector(0.0, 0.0, 1.0)
# Define the amplitude of the incident electric field
amplitude = norm(dir) # Amplitude of the incident electric field
# Create the incident field (uniform electric field)
ex = UniformField(; direction=dir , amplitude)
# Define the parameters
r_a = 0.0 #(r_in + r_PEC) / 2   # Radius of the PEC core
r_b = 1.0 #(r_Out + r_in) / 2  # Radius of the outer sphere

# Define the surrounding medium (outer PEC-Core )
ε_r1 = ε_r2 = rand()    # Relative permittivity of the surrounding medium (assuming free space) ε_r1 = ε_surrounding / ε0
ε_surrounding = ε_r1 * ε0  
medium_surrounding = Medium(ε_surrounding, μ0)

# Define the PEC material for the core (inner layer)
medium_pec_core = Medium(0.0 , μ0)  # PEC has zero permittivity and unit permeability

# define observation points in both layers and outside of the sphere
r_in = r_b / 2 
r_out =  r_b * 2

OB_Die = SVector(0.0, r_in, 0.0 )
OB_Out = SVector(r_out, 0.0, 0.0)

point_cart = [OB_Die , OB_Out] 



# Create the layered sphere with the PEC core
#sp = LayeredSphere(; radii=SVector(r_a, r_b), filling=SVector(medium_pec_core, medium_surrounding))
#sp = PECSphere(r_a)
sp = DielectricSphere(; radius= r_b, filling= medium_surrounding)
# Compute the scalar potential at the observation point
Φ = scatteredfield(sp, ex, ScalarPotential(point_cart))
Amplitude = 1.0
# the coefficient for inside and outside of Sphere with PEC_Core
l_max = 100
A_SC = calculate_LH_A2(Amplitude, r_a, r_b, ε_r1, ε_r2, l_max) # the coefficient for inside the SphrScatLib  or A_SC = calculate_ASph_with_Core(Amplitude, r_a, r_b, ε_r1, l_max)
E_SC = calculate_OH_E(Amplitude, r_b, ε_r1, ε_r2, l_max) # the coefficient for outside the SphrScatLib  or E_SC =  calculate_ESph_with_Core(Amplitude, r_b, ε_r1, l_max)

potential_in = A_SC[2] * r_in
potential_out = (E_SC[2] * r_out^(-2)) # - (amplitude * r_out)   

@test Φ[1]   ≈ potential_out 
@test Φ[2]   ≈ potential_out 


