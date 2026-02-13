module AtomicSplines

using LinearAlgebra
using SparseArrays
using FastGaussQuadrature
using Arpack

# Exportaciones
# export BSplineBasis, generate_basis
# export bspline, d_bspline, eval_expansion
# export assemble_core, assemble_J_matrix, assemble_centrifugal, assemble_operator_matrix, assemble_interaction_matrix, assemble_exchange_matrix, assemble_projection_operator
# export solve_poisson_potential, solve_poisson_general
# export solve_eigen, Orbital, Atom
# export solve_orbital!
# export wigner_3j_squared, generate_interactions, assemble_fock_matrix
# export Shell, AtomicConfig, InteractionTerm

# Basis
export BSplineBasis, eval_bspline_kernel!,  generate_basis
# Integrals
export assemble_core, assemble_J_matrix
# Poisson
export solve_poisson_potential
# Helpers
export evaluate_orbital


# Structs
include("structs.jl") 
include("config.jl")

# Math logic
include("basis.jl")

# Helpers
include("helpers.jl")

# Assemblers
include("integrals.jl")
include("poisson.jl")

# Linear algebra
include("solver.jl")

# High level logic
include("atomic_calculations.jl")


end
