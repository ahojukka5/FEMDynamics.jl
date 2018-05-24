# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/LinearImplicitDynamics/blob/master/LICENSE

using JuliaFEM
using JuliaFEM.Preprocess
using Logging
Logging.configure(level=DEBUG)
add_elements! = JuliaFEM.add_elements! # FIXME

# to load matrices:
# using JLD
# data = load("data.jld")
# M = data["mass matrix"]
# K = data["stiffness matrix"]
# mesh = data["mesh"]

# read mesh
mesh = abaqus_read_mesh("model.inp")
info("element sets = ", collect(keys(mesh.element_sets)))
info("surface sets = ", collect(keys(mesh.surface_sets)))
 
# create a field problem
beam = Problem(Elasticity, "BEAM", 3)
beam_elements = create_elements(mesh, "BEAM")
update!(beam_elements, "youngs modulus", 165.0)
update!(beam_elements, "poissons ratio", 0.275)
update!(beam_elements, "density", 7.10E-9)
add_elements!(beam, beam_elements) 

assemble!(beam, 0.0)
assemble_mass_matrix!(beam, 0.0)

# store matrices to h5 file format
using JLD
K = sparse(beam.assembly.K)
M = sparse(beam.assembly.M)
save("matrices.jld", "stiffness matrix", K, "mass matrix", M, "mesh", mesh; compress=true)
