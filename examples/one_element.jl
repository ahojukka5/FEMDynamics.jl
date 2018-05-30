# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/LinearImplicitDynamics.jl/blob/master/LICENSE

using JuliaFEM

X = Dict(1 => [0.0, 0.0],
         2 => [2.0, 0.0],
         3 => [2.0, 2.0],
         4 => [0.0, 2.0])

element = Element(Quad4, [1, 2, 3, 4])
element.id = 1
update!(element, "geometry", X)
update!(element, "youngs modulus", 288.0)
update!(element, "poissons ratio", 1/3)
update!(element, "density", 36.0)
problem = Problem(Elasticity, "2x2 element", 2)
problem.properties.formulation = :plane_stress
add_elements!(problem, [element])
time = 0.0
assemble!(problem, time)
assemble_mass_matrix!(problem, time)

K = sparse(problem.assembly.K)
M = sparse(problem.assembly.M)
cfM = cholfact(M)

function model1(dx, x, p, t)
    f = zeros(8)
    f[[1, 3, 5, 7]] = 5.0
    f[[2, 4, 6, 8]] = -10.0
    f[1] += 10
    ndofs = round(Int, length(dx)/2)
    u = x[1:ndofs]
    v = x[ndofs+1:end]
    dx[1:ndofs] = x[ndofs+1:end]
    dx[ndofs+1:end] = cfM \ (f - K*u)
end

using DifferentialEquations

ndofs = size(K, 1)
u0 = zeros(ndofs)
v0 = zeros(ndofs)
x0 = [u0; v0]
tspan = (0.0, 10.0)
prob = ODEProblem(model1, x0, tspan)
sol = solve(prob; dtmax=0.1)
u = hcat([x[1:ndofs] for x in sol.u]...)
v = hcat([x[ndofs+1:end] for x in sol.u]...)

xmid = mean(u[1:2:end, :], 1)
ymid = mean(u[2:2:end, :], 1)
vxmid = mean(v[1:2:end, :], 1)
vymid = mean(v[2:2:end, :], 1)

t0, t1 = tspan
nsteps = length(sol.u)
println("Number of time steps = $nsteps")
t = linspace(t0, t1, nsteps)

dim = 2
nnodes = length(X)
# ndofs = dim * nnodes

for (time, x) in zip(t, sol.u)
    u = x[1:ndofs]
    v = x[ndofs+1:end]
    ur = reshape(u, dim, nnodes)
    vr = reshape(v, dim, nnodes)
    ud = Dict(j => ur[:,j] for j=1:4)
    vd = Dict(j => vr[:,j] for j=1:4)
    update!(element, "displacement", time => ud)
    update!(element, "velocity", time => vd)
end

step = Analysis(Nonlinear)
add_problems!(step, [problem])

xdmf = Xdmf("one_element_results"; overwrite=true)
add_results_writer!(step, xdmf)

for time in t
    JuliaFEM.write_results!(step, time)
end

close(xdmf.hdf)
