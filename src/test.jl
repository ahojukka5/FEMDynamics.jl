# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/LinearImplicitDynamics.jl/blob/master/LICENSE

using DifferentialEquations

# As an example, I took this:

# http://docs.juliadiffeq.org/latest/tutorials/ode_example.html#Example-2:-Solving-Systems-of-Equations-1

function lorenz(du,u,p,t)
 du[1] = 10.0*(u[2]-u[1])
 du[2] = u[1]*(28.0-u[3]) - u[2]
 du[3] = u[1]*u[2] - (8/3)*u[3]
end

u0 = [1.0;0.0;0.0]
tspan = (0.0,100.0)
prob = ODEProblem(lorenz,u0,tspan)
sol = solve(prob)

# Our example Mu'' + Ku = f(t), dim(u) = 2

# 1. Build model matrices

ke = [1.0 -1.0; -1.0 1.0]
me = [2.0 1.0; 1.0 2.0]
f(t) = [0.0, 1.0]*sin.(t)

nelem = 2
ndofs = 2*nelem-1
K = zeros(ndofs, ndofs)
M = zeros(ndofs, ndofs)
for j=1:ndofs-1
    K[j:j+1, j:j+1] += ke
    M[j:j+1, j:j+1] += me
end

free_dofs = collect(2:ndofs)
K_red = K[free_dofs, free_dofs]
M_red = M[free_dofs, free_dofs]
invMK = M_red\K_red

# 2. Construct model

function model(dx,x,p,t)
    du = dx[1:2]
    dv = dx[3:4]
    u = x[1:2]
    v = x[3:4]
    x[1:2] = v
    x[3:4] = invMK*u + f(t)
end

# 3. Solve

x0 = zeros(4)
tspan = (0.0, 10.0)
prob = ODEProblem(model, x0, tspan)
sol = solve(prob; dtmax=0.1)

println("Ready. Lenght of solution vector: $(length(sol.u))")

# 4. Visualize solution

using Plots

u1 = [ui[1] for ui in sol.u]
u2 = [ui[2] for ui in sol.u]
t = linspace(tspan[1], tspan[2], length(sol.u))
plot(t, [u1, u2])
gui()

# Results are wrong.
