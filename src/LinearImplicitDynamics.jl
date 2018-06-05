# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/LinearImplicitDynamics.jl/blob/master/LICENSE

""" Linear, implicit dynamics solver for JuliaFEM. """
module LinearImplicitDynamics

using FEMBase
using DifferentialEquations

type LinearImplicit <: AbstractAnalysis
    tspan :: Tuple{Float64,Float64}
    sol
end

function LinearImplicit()
    tspan = (0.0, 1.0)
    return LinearImplicit(tspan, nothing)
end

function FEMBase.run!(analysis::Analysis{LinearImplicit})

    time = analysis.properties.tspan[1]
    assemble!(analysis, time; with_mass_matrix=true)

    M = SparseMatrixCOO()
    K = SparseMatrixCOO()
    Kg = SparseMatrixCOO()
    f = SparseMatrixCOO()
    fg = SparseMatrixCOO()

    for problem in get_problems(analysis)
        is_field_problem(problem) || continue
        append!(M, problem.assembly.M)
        append!(K, problem.assembly.K)
        append!(Kg, problem.assembly.Kg)
        append!(f, problem.assembly.f)
        append!(fg, problem.assembly.fg)
    end

    dim = size(K, 1)

    M = sparse(M, dim, dim)
    K = sparse(K, dim, dim) + sparse(Kg, dim, dim)
    f = sparse(f, dim, 1) + sparse(fg, dim, 1)

    for problem in get_problems(analysis)
        is_boundary_problem(problem) || continue
        eliminate_boundary_conditions!(problem, K, M, f)
    end

    K = 1/2*(K + K')
    M = 1/2*(M + M')
    SparseArrays.droptol!(K, 1.0e-9)
    SparseArrays.droptol!(M, 1.0e-9)
    nz = get_nonzero_rows(K)
    K_red = K[nz,nz]
    M_red = M[nz,nz]
    cfM = cholfact(M_red)

    function model(dx, x, p, t)
        ndofs = round(Int, length(dx)/2)
        u = x[1:ndofs]
        ux = round(mean(u[1:3:end]), 2)
        uy = round(mean(u[2:3:end]), 2)
        uz = round(mean(u[3:3:end]), 2)
        println("$t: coords = ($ux, $uy, $uz)")
        v = x[ndofs+1:end]
        dx[1:ndofs] = x[ndofs+1:end]
        dx[ndofs+1:end] = cfM \ (f - K_red*u)
    end

    ndofs = size(K_red, 1)
    u0 = zeros(ndofs)
    v0 = zeros(ndofs)
    v0[1:3:end] += 20.0
    v0[2:3:end] += 10.0
    x0 = [u0; v0]
    tspan = analysis.properties.tspan
    prob = ODEProblem(model, x0, tspan)
    analysis.properties.sol = solve(prob; dtmax=0.1)
    #u = hcat([x[1:ndofs] for x in sol.u]...)
    #v = hcat([x[ndofs+1:end] for x in sol.u]...)
end

export LinearImplicit

end
