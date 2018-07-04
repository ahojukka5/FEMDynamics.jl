# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMDynamics.jl/blob/master/LICENSE

using FEMDynamics
using JuliaFEM
using FEMBase.Test

X = Dict(1 => [0.0, 0.0],
         2 => [2.0, 0.0],
         3 => [2.0, 2.0],
         4 => [0.0, 2.0])

element = Element(Quad4, [1, 2, 3, 4])
update!(element, "geometry", X)
update!(element, "youngs modulus", 288.0)
update!(element, "poissons ratio", 1/3)
update!(element, "density", 36.0)
update!(element, "displacement load 2", 36.0*9.81)
problem = Problem(Elasticity, "2x2 element", 2)
problem.properties.formulation = :plane_stress
add_elements!(problem, [element])
analysis = Analysis(LinearImplicit, "dropping element")
add_problems!(analysis, [problem])
run!(analysis)
midpnt_u = element("displacement", (0.0, 0.0), Inf)
midpnt_v = element("velocity", (0.0, 0.0), Inf)
@test isapprox(midpnt_u, [0.0, 1/2*9.81])
@test isapprox(midpnt_v, sqrt.(2.0*9.81*midpnt_u))
