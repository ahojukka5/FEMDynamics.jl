# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/LinearImplicitDynamics.jl/blob/master/LICENSE

using FEMBase
using LinearImplicitDynamics

using Base.Test

@testset "LinearImplicitDynamics.jl" begin
    # Test stiffness matrix
    K = 1.0
    K_expected = 1.0
    @test isapprox(K, K_expected)
end
