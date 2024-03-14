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
r_a = 0.0 #(r_in + r_PEC) / 2   # Radius of the PEC core
r_b = 1.0 #(r_Out + r_in) / 2  # Radius of the outer sphere

# Define the surrounding medium (outer PEC-Core )
ε_r1 = ε_r2 = 5 #rand()    # Relative permittivity of the surrounding medium (assuming free space) ε_r1 = ε_surrounding / ε0
ε_surrounding = ε_r1 * ε0  
medium_surrounding = Medium(ε_surrounding, μ0)

# Define the PEC material for the core (inner layer)
medium_pec_core = Medium(0.0 , μ0)  # PEC has zero permittivity and unit permeability

# define observation points in both layers and outside of the sphere
r_in = r_b / 2 
r_out =  r_b * 2

θ = rand()

OB_Die = SVector(0.0, r_in * sin(θ), r_in * cos(θ) )
OB_Out = SVector(0.0, r_out * sin(θ), r_out * cos(θ))

point_cart = [OB_Die , OB_Out] 



# Create the layered sphere with the PEC core
#sp = LayeredSphere(; radii=SVector(r_a, r_b), filling=SVector(medium_pec_core, medium_surrounding))
#sp = PECSphere(r_a)
sp = DielectricSphere(; radius= r_b, filling= medium_surrounding)

# the coefficient for inside and outside of Sphere with PEC_Core
l_max = 100
A_SC = calculate_LH_A2(amplitude, r_a, r_b, ε_r1, ε_r2, l_max) # the coefficient for inside the SphrScatLib  or A_SC = calculate_ASph_with_Core(Amplitude, r_a, r_b, ε_r1, l_max)
E_SC = calculate_OH_E(amplitude, r_b, ε_r1, ε_r2, l_max) # the coefficient for outside the SphrScatLib  or E_SC =  calculate_ESph_with_Core(Amplitude, r_b, ε_r1, l_max)


function calculate_Phi_2(amplitude, r, ra, rb, xi, l_max)
    Phi_2 = 0.0
    Al = calculate_LH_A2(amplitude, ra, rb, ε_r1, ε_r2, l_max)
    for l in 0:l_max-1
        term = Al[l+1] * (r^l - ra^(2l+1) * r^-(l+1)) * Pl(xi,l)
        Phi_2 += term
    end
    return Phi_2
end

# Compute the scalar potential at the observation point
Φ = scatteredfield(sp, ex, ScalarPotential(point_cart))
potential_in = calculate_Phi_2(amplitude, r_in, r_a, r_b, cos(θ) , l_max) 
potential_out = ((E_SC[2] * r_out^(-2))  - (amplitude * r_out)) * cos(θ)  
@test Φ[1]  - ((amplitude * r_in) * cos(θ))  ≈ potential_in
@test Φ[2]  - ((amplitude * r_out) * cos(θ))  ≈ potential_out 


