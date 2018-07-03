# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/LinearImplicitDynamics.jl/blob/master/LICENSE

""" Linear, implicit dynamics solver for JuliaFEM. """
module LinearImplicitDynamics

using Reexport
@reexport using FEMBase
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

    for problem in get_problems(analysis)
        assemble!(problem, time)
        assemble_mass_matrix!(problem, time)
    end

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
    dropzeros!(K)
    dropzeros!(M)
    nz = get_nonzero_rows(K)
    K_red = K[nz,nz]
    M_red = M[nz,nz]
    cfM = cholfact(M_red)
    display(full(K_red))
    display(full(M_red))

    function model(dx, x, p, t)
        ndofs = round(Int, length(dx)/2)
        u = x[1:ndofs]
        println("Solving at time $t")
        v = x[ndofs+1:end]
        dx[1:ndofs] = x[ndofs+1:end]
        dx[ndofs+1:end] = cfM \ (f - K_red*u)
    end

    ndofs = size(K_red, 1)
    u0 = zeros(ndofs)
    v0 = zeros(ndofs)
    x0 = [u0; v0]
    prob = ODEProblem(model, x0, analysis.properties.tspan)
    sol = analysis.properties.sol = solve(prob)

    # FIXME: support for variable number of dofs / node
    dim = get_unknown_field_dimension(first(get_problems(analysis)))
    nnodes = Int(ndofs/dim)

    info("Number of time steps = ", length(sol))

    function to_dict(u, dim, nnodes)
        return Dict(j => [u[dim*(j-1)+k] for k=1:dim] for j=1:nnodes)
    end

    for (time, x) in zip(sol.t, sol.u)
        u = to_dict(x[1:ndofs], dim, nnodes)
        v = to_dict(x[ndofs+1:end], dim, nnodes)
        for problem in get_problems(analysis)
            for elements in get_elements(problem)
                update!(elements, "displacement", time => u)
                update!(elements, "velocity", time => v)
            end
        end
    end

end

# function FEMBase.write_results!(analysis::LinearImplicitDynamics, writer::Xdmf)
#     # Not implemented yet
# end

export LinearImplicit

end
