# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMDynamics.jl/blob/master/LICENSE

""" Dynamics solver for JuliaFEM. """
module FEMDynamics

using Reexport
@reexport using FEMBase
using DifferentialEquations

mutable struct LinearImplicit <: AbstractAnalysis
    tspan :: Tuple{Float64,Float64}
    sol
    K :: SparseMatrixCSC{Float64}
    M :: SparseMatrixCSC{Float64}
    f :: SparseVector{Float64}
    u0 :: SparseVector{Float64}
    v0 :: SparseVector{Float64}
end

function LinearImplicit()
    tspan = (0.0, 1.0)
    return LinearImplicit(tspan, nothing, spzeros(0,0), spzeros(0,0),
                          spzeros(0), spzeros(0), spzeros(0))
end

"""
    get_global_matrices(analysis)

Assemble global matrices for problem.
"""
function get_global_matrices(analysis::Analysis{LinearImplicit})
    p = analysis.properties
    if !isempty(p.K) && !isempty(p.M)
        return p.K, p.M, p.f
    end

    time = first(analysis.properties.tspan)

    for problem in get_problems(analysis)
        assemble!(problem, time)
        is_field_problem(problem) && assemble_mass_matrix!(problem, time)
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

    p.K = 1/2*(K + K')
    p.M = 1/2*(M + M')
    p.f = f
    return p.K, p.M, p.f
end

function FEMBase.run!(analysis::Analysis{LinearImplicit})

    # Preprocess: create global matrices, eliminate boundary conditions, remove
    # any nonzero rows and factorize

    K, M, f = get_global_matrices(analysis)
    ndofs = size(K, 1)
    full_vec = zeros(ndofs) # for temporary operations
    dim = get_unknown_field_dimension(first(get_problems(analysis)))
    nnodes = Int(ndofs/dim)

    for problem in get_problems(analysis)
        is_boundary_problem(problem) || continue
        eliminate_boundary_conditions!(problem, K, M, f)
    end
    dropzeros!(K)
    dropzeros!(M)

    if isempty(analysis.properties.u0)
        analysis.properties.u0 = spzeros(ndofs)
    end
    if isempty(analysis.properties.v0)
        analysis.properties.v0 = spzeros(ndofs)
    end
    u0 = analysis.properties.u0
    v0 = analysis.properties.v0

    nz = get_nonzero_rows(K)
    K = K[nz,nz]
    M = M[nz,nz]
    f = f[nz]
    u0 = u0[nz]
    v0 = v0[nz]
    M_fact = cholfact(M)
    ndofs = size(K, 1)

    # Solution

    function model(dx, x, p, t)
        info("Solving at time $t")
        ndofs = round(Int, length(dx)/2)
        u = x[1:ndofs]
        v = x[ndofs+1:end]
        dx[1:ndofs] = v
        dx[ndofs+1:end] = M_fact \ full(f - K*u)
    end

    x0 = [u0; v0]
    prob = ODEProblem(model, x0, analysis.properties.tspan)
    sol = analysis.properties.sol = solve(prob)

    # postprocess, update displacement and velocity back to elements

    # FIXME: support for variable number of dofs / node
    function to_dict(u, dim, nnodes)
        return Dict(j => [u[dim*(j-1)+k] for k=1:dim] for j=1:nnodes)
    end

    info("Number of nodes = $nnodes, number of time steps = $(length(sol))")

    for (time, x) in zip(sol.t, sol.u)
        full_vec[nz] = x[1:ndofs]
        u = to_dict(full_vec, dim, nnodes)
        full_vec[nz] = x[ndofs+1:end]
        v = to_dict(full_vec, dim, nnodes)
        for problem in get_problems(analysis)
            is_field_problem(problem) || continue
            for elements in get_elements(problem)
                update!(elements, "displacement", time => u)
                update!(elements, "velocity", time => v)
            end
        end
    end

end

# function FEMBase.write_results!(analysis::ImplicitDynamics, writer::Xdmf)
#     # Not implemented yet
# end

export LinearImplicit

end
