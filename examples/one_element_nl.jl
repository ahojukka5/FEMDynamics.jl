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
try
    xdmf = Xdmf("one_element_results_nl"; overwrite=true)
catch
    close(xdmf.hdf)
    xdmf = Xdmf("one_element_results_nl"; overwrite=true)
end
#add_results_writer!(step, xdmf)

time = 0.0
dt = 0.001
sr = 0.9 # spectral radius
a_m = (2*sr-1)/(sr+1)
a_f = sr/(sr+1)
beta = 1/4*(1 - a_m + a_f)^2
gamma = 1/2 - a_m + a_f
nnodes = length(X)
dim = 2
ndofs = nnodes * dim
d = Dict(1 => zeros(ndofs))
v = Dict(1 => zeros(ndofs))
#v[1][[2*(2-1)+2, 2*(3-1)+2]] = -150.0
#v[1][[2*(5-1)+2, 2*(6-1)+2]] = -100.0
a = Dict(1 => zeros(ndofs))
E_int = Dict(1 => 0.0)
E_kin = Dict(1 => 0.0)
E_tot = Dict(1 => 0.0)
du = zeros(ndofs)

assemble_mass_matrix!(problem, 0.0)

function to_dict(u, dim, nnodes)
    return Dict(j => [u[dim*(j-1)+k] for k=1:dim] for j=1:nnodes)
end

update!(elements, "displacement", 0.0 => to_dict(d[1], dim, nnodes))
JuliaFEM.write_results!(step, time)

for n=1:600
    time = n*dt
    # info("Starting time step $n: $time")
    d[n+1] = copy(d[n])
    for i=1:10
        # info("Starting nonlinear iteration $i")
        empty!(problem.assembly)
        assemble!(problem, time)
        K = full(problem.assembly.K)
        M = full(problem.assembly.M)
        r = -full(problem.assembly.f)
        C = zeros(K)

        v[n+1] = gamma/(beta*dt)*(d[n+1]-d[n]) - (gamma-beta)/beta*v[n] - (gamma-2.0*beta)/(2.0*beta)*dt*a[n]
        a[n+1] = 1.0/(beta*dt^2)*(d[n+1]-d[n]) - 1.0/(beta*dt)*v[n] - (1.0-2.0*beta)/(2.0*beta)*dt*a[n]

        # generalized alpha
        d_a = (1-a_f)*d[n+1] + a_f*d[n]
        v_a = (1-a_f)*v[n+1] + a_f*v[n]
        a_a = (1-a_m)*a[n+1] + a_m*a[n]
        r_a = (1-a_f)*r + a_f*r

        K_effdyn = (1.0-a_m)/(beta*dt^2)*M + (1.0-gamma)/(beta*dt)*C + (1-a_f)*K
        r_effdyn = M*a_a + C*v_a + r_a

        fill!(du, 0.0)
        fixed_dofs = [1, 2]
        free_dofs = setdiff(collect(1:ndofs), fixed_dofs)
        du[free_dofs] = K_effdyn[free_dofs, free_dofs] \ -r_effdyn[free_dofs]
        d[n+1] = d[n+1] + du
        update!(elements, "displacement", time => to_dict(d[n+1],dim,nnodes))
        update!(elements, "velocity", time => to_dict(v[n+1],dim,nnodes))
        update!(elements, "acceleration", time => to_dict(a[n+1],dim,nnodes))
        E_kin[n+1] = 1/2*dot(v[n+1], M*v[n+1])
        E_int[n+1] = 1/2*dot(d[n+1], K*d[n+1])
        E_tot[n+1] = 1/2*dot(d[n+1], K_effdyn*d[n+1]) - dot(r_effdyn, d[n+1])
        dunorm = norm(du)
        if dunorm < 1.0e-6
            info("n=$n: nonlinear iteration converged in $i iterations.")
            break
        end
    end
    #JuliaFEM.write_results!(step, time)
end

#close(xdmf.hdf)

using Plots
s = sort(collect(keys(E_int)))
#E_tot = Dict(j => E_int[j]+E_kin[j] for j in s)
E1 = [E_int[j] for j in s]
E2 = [E_kin[j] for j in s]
E3 = [E_tot[j] for j in s]
plt = plot(s, E1, label="internal energy", c=:red, w=2)
plot!(plt, s, E2, label="kinetic energy", c=:blue, w=2)
plot!(plt, s, E3, label="total energy", c=:green, w=2)
display(plt)
