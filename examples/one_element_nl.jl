# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/LinearImplicitDynamics.jl/blob/master/LICENSE

using JuliaFEM
using Logging
Logging.configure(level=DEBUG)

X = Dict(1 => [0.0, 0.0],
         2 => [5.0, 0.0],
         3 => [4.0, 1.0],
         4 => [0.0, 1.0],
         5 => [5.0, 5.0],
         6 => [4.0, 5.0])

element1 = Element(Quad4, [1, 2, 3, 4])
element2 = Element(Quad4, [2, 5, 6, 3])
elements = [element1, element2]
update!(elements, "geometry", X)
update!(elements, "youngs modulus", 5.0)
update!(elements, "poissons ratio", 0.3)
update!(elements, "density", 1.0e-5)
update!(elements, "displacement load 2", -0.1)
problem = Problem(Elasticity, "2x2 element", 2)
problem.properties.formulation = :plane_stress
problem.properties.finite_strain = true
problem.properties.geometric_stiffness = true
add_elements!(problem, elements)
step = Analysis(Nonlinear)
add_problems!(step, [problem])
xdmf = Xdmf("one_element_results_nl"; overwrite=true)
add_results_writer!(step, xdmf)

time = 0.0
dt = 0.001
beta = 1/4
gamma = 1/2
nnodes = length(X)
dim = 2
ndofs = nnodes * dim
d = Dict(1 => zeros(ndofs))
v = Dict(1 => zeros(ndofs))
a = Dict(1 => zeros(ndofs))
du = zeros(ndofs)

assemble_mass_matrix!(problem, 0.0)

function to_dict(u, dim, nnodes)
    return Dict(j => [u[dim*(j-1)+k] for k=1:dim] for j=1:nnodes)
end

update!(elements, "displacement", 0.0 => to_dict(d[1], dim, nnodes))
JuliaFEM.write_results!(step, time)

for n=1:600
    time = n*dt
    info("Starting time step $n: $time")
    d[n+1] = copy(d[n])
    for i=1:10
        info("Starting nonlinear iteration $i")
        empty!(problem.assembly)
        assemble!(problem, time)
        K = full(problem.assembly.K)
        M = full(problem.assembly.M)
        r = -full(problem.assembly.f)
        C = zeros(K)
        K_effdyn = 1.0/(beta*dt^2)*M + gamma/(beta*dt)*C + K
        v[n+1] = gamma/(beta*dt)*(d[n+1]-d[n]) - (gamma-beta)/beta*v[n] - (gamma-2.0*beta)/(2.0*beta)*dt*a[n]
        a[n+1] = 1.0/(beta*dt^2)*(d[n+1]-d[n]) - 1.0/(beta*dt)*v[n] - (1.0-2.0*beta)/(2.0*beta)*dt*a[n]
        r_effdyn = M*a[n+1] + C*v[n+1] + r

        fill!(du, 0.0)
        fixed_dofs = [1, 2]
        free_dofs = setdiff(collect(1:ndofs), fixed_dofs)
        du[free_dofs] = K_effdyn[free_dofs, free_dofs] \ -r_effdyn[free_dofs]
        d[n+1] = d[n+1] + du
        update!(elements, "displacement", time => to_dict(d[n+1],dim,nnodes))
        update!(elements, "velocity", time => to_dict(v[n+1],dim,nnodes))
        update!(elements, "acceleration", time => to_dict(a[n+1],dim,nnodes))
        dunorm = norm(du)
        info("||du|| = $dunorm")
        if dunorm < 1.0e-6
            info("Nonlinear iteration converged.")
            break
        end
    end
    JuliaFEM.write_results!(step, time)
end

close(xdmf.hdf)
