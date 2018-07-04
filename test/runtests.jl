# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMDynamics.jl/blob/master/LICENSE

using FEMDynamics
using FEMBase.Test

@testset "LinearImplicit" begin
    include("test_dropping_element.jl")
end
