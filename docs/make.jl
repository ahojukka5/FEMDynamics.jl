# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/LinearImplicitDynamics.jl/blob/master/LICENSE

using Documenter, LinearImplicitDynamics

makedocs(modules=[LinearImplicitDynamics],
         format = :html,
         checkdocs = :all,
         sitename = "LinearImplicitDynamics.jl",
         pages = ["index.md"]
        )
