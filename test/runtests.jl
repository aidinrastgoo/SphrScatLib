using Revise
using Test
using SphrScatLib

@testset "SphrScatLib" begin
    
    include("compare_Sphr.jl")
    include("compare_Sphr_Core.jl")
    include("SphericalScattering.jl")

end