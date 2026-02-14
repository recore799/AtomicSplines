module AtomicSplines

using LinearAlgebra
using SparseArrays
using FastGaussQuadrature
using Arpack

# Basis
export BSplineBasis, SolverWorkspace, generate_basis, eval_bspline_kernel!
# Integrals
export assemble_core, assemble_J_matrix, init_scf_workspace, build_interaction_tensors, assemble_J_matrix_param, assemble_J_matrix_tensor
# Poisson
export solve_poisson_potential, solve_poisson_fast, solve_poisson_tensor
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
