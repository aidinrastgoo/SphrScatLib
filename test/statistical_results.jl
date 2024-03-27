using Plots
using SpecialFunctions  # For Legendre polynomials
using LegendrePolynomials
ε_r1 = rand(1:99)
ε_r2 = rand(1:99)
ξ= rand()
rb = 1.0 # radius_b take as refrence radius and has value of 1.0 p.u 
ra = 0.5 # radius_a is half of radius_a
radius = (50 + rand(1:50)) / 100
l_max = 30 
amplitude =1.0# 


function calculate_Phi_1(amplitude, r, ra, rb, ξ, l_max)
    Al = calculate_HH_A1(amplitude, ra, rb, ε_r1, ε_r2, l_max)
    Phi_1 = zeros(l_max)
    Phi_1[1] = Al[1]
    for l in 1:l_max-1
        Phi_1[l+1] = - Al[l+1] * ((r^l) - ((r^(-(l+1)))* ra^(2l+1)))  * Pl(ξ,l) 
    end
    return sum(Phi_1)
end

# Calculate Potential of northern hiemspere for each value of l
φ_values = zeros(l_max)
for l in 2:l_max
    φ = calculate_Phi_1(amplitude, radius , ra, rb, ξ , l)
    φ_values[l] = φ  # Stor elements
end

# Plot the results
p = plot(1:l_max, φ_values , xlabel="l", ylabel="φ", label="", title="φ2 Vs. l")
annotate!([(l_max, maximum(φ_values)-(60*maximum(φ_values)/100), text("ε_r1 = $ε_r1", :right, 8)),
            (l_max, maximum(φ_values)-(70*maximum(φ_values)/100), text("ε_r2 = $ε_r2", :right, 8)),
            (l_max, maximum(φ_values)-(80*maximum(φ_values)/100), text("Radius = $radius", :right, 8)),
            (l_max, maximum(φ_values)-(90*maximum(φ_values)/100), text("cos(θ) = $ξ", :right, 8))])

# Display the  plot
display(p)