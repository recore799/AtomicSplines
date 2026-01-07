module AtomicSplines

using LinearAlgebra
using SparseArrays
using FastGaussQuadrature

export BSplineBasis, generate_basis
export bspline, d_bspline, eval_expansion
export assemble_core, assemble_J_matrix
export solve_poisson_potential

include("basis.jl")
include("grid.jl")
include("integrals.jl")
include("poisson.jl")

end
