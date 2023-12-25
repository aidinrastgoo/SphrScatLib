include("compute_coefficient_lower_hemisphere.jl")
include("compute_η.jl")

function calculate_HH_A1(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)
    A2 = calculate_LH_A2(Amplitude, radius_a, radius_b, ε_r1, ε_r2, l_max)
    η= calculate_η(l_max, ε_r1, ε_r2)
    return η.* A2
end