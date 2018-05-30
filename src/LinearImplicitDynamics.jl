# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/LinearImplicitDynamics.jl/blob/master/LICENSE

""" Linear, implicit dynamics solver for JuliaFEM. """
module LinearImplicitDynamics

using FEMBase

import FEMBase: solve!

type LinearImplicit <: AbstractSolver
end

function run!(solver::Solver{LinearImplicit}, time)
    # run analysis
end

export LinearImplicit

end
