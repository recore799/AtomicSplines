module AtomicSplines

using LinearAlgebra
using SparseArrays
using FastGaussQuadrature
using Arpack

# Exportaciones
export BSplineBasis, generate_basis
export bspline, d_bspline, eval_expansion
export assemble_core, assemble_J_matrix, assemble_centrifugal
export solve_poisson_potential
export solve_eigen, Orbital, Atom
export solve_orbital!

# Structs
include("structs.jl") 

# Math logic
include("basis.jl")
include("grid.jl")

# Assemblers
include("integrals.jl")
include("poisson.jl")

# Linear algebra
include("solver.jl")

# High level logic
include("atomic_calculations.jl")

end
