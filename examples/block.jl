# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/LinearImplicitDynamics.jl/blob/master/LICENSE

using JuliaFEM
using JuliaFEM.Preprocess
add_elements! = JuliaFEM.add_elements!
using LinearImplicitDynamics
using Logging
Logging.configure(level=DEBUG)

## Generating structured block mesh

# TODO: this should be in some package!
function quad_mesh(nel_x=10, nel_y=10; lx=1.0, ly=1.0)
    nnodes_x = nel_x+1
    nnodes_y = nel_y+1
    nnode = nnodes_x*nnodes_y
    nodemap = reshape(1:nnode, nnodes_x, nnodes_y)
    nodes_1 = vec(nodemap[1:nnodes_x-1, 1:nnodes_y-1])
    nodes_2 = vec(nodemap[2:nnodes_x, 1:nnodes_y-1])
    nodes_3 = vec(nodemap[2:nnodes_x, 2:nnodes_y])
    nodes_4 = vec(nodemap[1:nnodes_x-1, 2:nnodes_y])

    mesh = Mesh()

    nid = 1
    for y in linspace(0, ly, nnodes_y)
        for x in linspace(0, lx, nnodes_x)
            add_node!(mesh, nid, [x, y])
            add_node_to_node_set!(mesh, :NALL, nid)
            nid += 1
        end
    end

    elid = 1
    for c in zip(nodes_1, nodes_2, nodes_3, nodes_4)
        add_element!(mesh, elid, :Quad4, collect(c))
        add_element_to_element_set!(mesh, :EALL, elid)
        elid += 1
    end

    return mesh
end

nel_x = 10
nel_y = 3
mesh = quad_mesh(nel_x, nel_y; lx=5000.0, ly=90.0)
for nid=1:nel_x+1:(nel_x+1)*(nel_y+1)
    add_node_to_node_set!(mesh, :LEFT, nid)
end

block_elements = create_elements(mesh, "EALL")
update!(block_elements, "youngs modulus", 210.0e3)
update!(block_elements, "poissons ratio", 0.3)
update!(block_elements, "density", 7.80e-7)
update!(block_elements, "displacement load 2", 7.80e-7*9.81)
block = Problem(Elasticity, "10x1 block", 2)
block.properties.formulation = :plane_stress
add_elements!(block, block_elements)

fixed = Problem(Dirichlet, "fixed", 2, "displacement")
fixed_nodes = mesh.node_sets[:LEFT]
fixed_elements = [Element(Poi1, [j]) for j in fixed_nodes]
update!(fixed_elements, "geometry", mesh.nodes)
for j=1:2
    update!(fixed_elements, "displacement $j", 0.0)
end
add_elements!(fixed, fixed_elements)

xdmf = Xdmf("block_resulst"; overwrite=true)

# analysis = Analysis(Nonlinear, "modal analysis of block under self-weight")
# add_results_writer!(analysis, xdmf)
# add_problems!(analysis, [block, fixed])
# run!(analysis)

# analysis = Analysis(Modal, "modal analysis of block under self-weight")
# add_results_writer!(analysis, xdmf)
# add_problems!(analysis, [block, fixed])
# run!(analysis)

analysis = Analysis(LinearImplicit, "transient analysis of block under self-weight")
add_results_writer!(analysis, xdmf)
add_problems!(analysis, [block, fixed])
run!(analysis)

close(xdmf.hdf)

max_el_id = maximum(keys(mesh.elements))
element = first(filter(e -> e.id == max_el_id, block_elements))
disp = element.fields["displacement"].data
t = [d.first for d in disp]
u = [d.second[end][2] for d in disp]

# using Plots
# plot(t[1:10:end], u[1:10:end])
