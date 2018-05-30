# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/LinearImplicitDynamics.jl/blob/master/LICENSE

using JuliaFEM
using JuliaFEM.Preprocess
using JuliaFEM.Postprocess
using Logging
Logging.configure(level=INFO)
add_elements! = JuliaFEM.add_elements!

datadir = Pkg.dir("LinearImplicitDynamics", "examples", "ball")
meshfile = joinpath(datadir, "ball.med")
mesh = aster_read_mesh(meshfile)
mesh.element_sets[:BALL] = bset = Set{Int64}()
for (elid, eltype) in mesh.element_types
    if eltype == :Tet4
        push!(bset, elid)
    end
end

for (elset_name, element_ids) in mesh.element_sets
    nel = length(element_ids)
    println("Element set $elset_name contains $nel elements.")
end
for (nset_name, node_ids) in mesh.node_sets
    nno = length(node_ids)
    println("Node set $nset_name contains $nno nodes.")
end
nnodes = length(mesh.nodes)
println("Total number of nodes in mesh: $nnodes")
nelements = length(mesh.elements)
println("Total number of elements in mesh: $nelements")

ball_elements = create_elements(mesh, "BALL")
nel_ball_elements = length(ball_elements)
println("ball contains $nel_ball_elements elements")
update!(ball_elements, "youngs modulus", 288.0)
update!(ball_elements, "poissons ratio", 1/3)
update!(ball_elements, "density", 36.0e-3)
#update!(ball_elements, "displacement load 1", 0.20)
update!(ball_elements, "displacement load 2", 0.20)

xmax = -Inf
nid = 0
for (j, coord) in mesh.nodes
    x, y, z = coord
    if x > xmax
        xmax = x
        nid = j
    end
end
println("max_x node = $nid, xval = $xmax")

load_element = Element(Poi1, [nid])
update!(load_element, "geometry", mesh.nodes)

V_ball = 0.0
m_ball = 0.0
time = 0.0
for element in ball_elements
    for ip in get_integration_points(element)
        detJ = element(ip, time, Val{:detJ})
        rho = element("density", ip, time)
        V_ball += ip.weight * detJ
        m_ball += ip.weight * rho * detJ
    end
end
println("Ball volume: $V_ball")
println("Ball mass: $m_ball")

ball = Problem(Elasticity, "ball", 3)
add_elements!(ball, ball_elements)
update!(load_element, "displacement traction force z", 10e4)
# add_elements!(ball, [load_element])
empty!(ball.assembly)
assemble!(ball, time)
assemble_mass_matrix!(ball, time)

K = sparse(ball.assembly.K)
M = sparse(ball.assembly.M)
f = sparse(ball.assembly.f)
cfM = cholfact(M)

function model1(dx, x, p, t)
    k = 1.0 - 2.0/5.0*t
    #dim = get_unknown_field_dimension(ball) # = 3
    ndofs = round(Int, length(dx)/2)
    u = x[1:ndofs]
    ux = round(mean(u[1:3:end]), 2)
    uy = round(mean(u[2:3:end]), 2)
    uz = round(mean(u[3:3:end]), 2)
    println("$t: coords = ($ux, $uy, $uz)")
    v = x[ndofs+1:end]
    dx[1:ndofs] = x[ndofs+1:end]
    f2 = copy(f)
    #f2[3:3:end] *= k
    #f2[3*(nid-1)+2] += 100.0
    dx[ndofs+1:end] = cfM \ (k*f2 - K*u)
end

using DifferentialEquations

ndofs = size(K, 1)
u0 = zeros(ndofs)
v0 = zeros(ndofs)
v0[1:3:end] += 20.0
v0[2:3:end] += 10.0
x0 = [u0; v0]
tspan = (0.0, 10.0)
prob = ODEProblem(model1, x0, tspan)
sol = solve(prob; dtmax=0.1)
u = hcat([x[1:ndofs] for x in sol.u]...)
v = hcat([x[ndofs+1:end] for x in sol.u]...)

t0, t1 = tspan
nsteps = length(sol.u)
println("Number of time steps = $nsteps")
t = linspace(t0, t1, nsteps)

dim = 3
nnodes = length(mesh.nodes)
# ndofs = dim * nnodes
#
for (time, x) in zip(t, sol.u)
    u = x[1:ndofs]
    v = x[ndofs+1:end]
    ur = reshape(u, dim, nnodes)
    vr = reshape(v, dim, nnodes)
    ud = Dict(j => ur[:,j] for j=1:nnodes)
    vd = Dict(j => vr[:,j] for j=1:nnodes)
    update!(ball_elements, "displacement", time => ud)
    update!(ball_elements, "velocity", time => vd)
end

nid2 = find_nearest_node(mesh, [0.0, 0.0, 0.0])

ball("displacement", 0.0)
trajectory = [ball("displacement", time)[nid2] for time in t]
traj_x = [t[1] for t in trajectory]
traj_y = [t[2] for t in trajectory]
traj_z = [t[3] for t in trajectory]
xmin, xmax = minimum(traj_x), maximum(traj_x)
ymin, ymax = minimum(traj_y), maximum(traj_y)
zmin, zmax = minimum(traj_z), maximum(traj_z)
println("bounding box x = ($xmin, $xmax)")
println("bounding box y = ($ymin, $ymax)")
println("bounding box z = ($zmin, $zmax)")

step = Analysis(Nonlinear)
add_problems!(step, [ball])

xdmf = Xdmf("ball_results"; overwrite=true)
add_results_writer!(step, xdmf)

for time in linspace(t0, t1, 600)
    JuliaFEM.write_results!(step, time)
end

close(xdmf.hdf)
