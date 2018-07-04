# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMDynamics.jl/blob/master/LICENSE

using Documenter, FEMDynamics

makedocs(modules=[FEMDynamics],
         format = :html,
         checkdocs = :all,
         sitename = "FEMDynamics.jl",
         pages = ["index.md"]
        )
