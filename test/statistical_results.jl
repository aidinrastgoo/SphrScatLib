using Plots
ε_r1 = rand(1:99)
ε_r2 = rand(1:99)
ξ= rand()
rb = 1.0 # radius_b take as refrence radius and has value of 1.0 p.u 
ra = 0.5 # radius_a is half of radius_a
radius = (50 + rand(10:50)) / 100
l_max = 10 
amplitude =1.0# 


# Calculate Potential of northern hiemspere for each value of l
φ_in = zeros(l_max)
for l in 2:l_max
    φ = calculate_Phi(amplitude, radius , ra, rb, ξ ,ε_r1 ,ε_r2 , l)
    φ_in[l] = - φ  # Stor elements
end

# Plot the results
p = plot(1:l_max, φ_in, xlabel="values of l", ylabel="φ", label="", title="Potential Comparison in Hemisphere for adjacent values of l")

# Define font size for the title
titlefontsize = 5

title!("Potential Comparison in Hemisphere for adjacent values of l", titlefont=font(titlefontsize))

annotate!([(l_max, maximum(φ_in)-(60*maximum(φ_in)/100), text("ε_r1 = $ε_r1", :right, 8)),
            (l_max, maximum(φ_in)-(70*maximum(φ_in)/100), text("ε_r2 = $ε_r2", :right, 8)),
            (l_max, maximum(φ_in)-(80*maximum(φ_in)/100), text("Radius = $radius", :right, 8)),
            (l_max, maximum(φ_in)-(90*maximum(φ_in)/100), text("cos(θ) = $ξ", :right, 8))])

# Display the  plot
display(p)



l_max = 11
Phi = zeros(Float64, l_max)
for l in 1 : l_max
    r = 0.0 + (l * 0.1)
    phi = calculate_Phi(1.0, r , 0.5, 1.0 , 1.0 , 2 , 5 , 80)
    Phi[l] =  phi  # Stor elements
end

c = plot(1:l_max, Phi, xlabel="length", ylabel="φ", label="Potential", title="Potential Variation", titlefont=font(titlefontsize))
scatter!([1:4], [Phi[1:4]], markershape=:square, markersize=4, markercolor=:red, label="PEC Region")
scatter!([5:9], [Phi[5:9]], markershape=:circle, markersize=4, markercolor=:blue, label="Region 1")
scatter!([10:11], [Phi[10:11]], markershape=:square, markersize=4, markercolor=:green, label="Enviroment")

display(plot(p, c))