# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/LinearImplicitDynamics.jl/blob/master/LICENSE

using LinearImplicitDynamics
using FEMBase.Test

@testset "LinearImplicitDynamics.jl" begin
    include("test_dropping_element.jl")
end
