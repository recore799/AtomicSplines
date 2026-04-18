module AtomicSplines

using LinearAlgebra
using BandedMatrices
using FastGaussQuadrature
using StaticArrays

# Basis
export BSplineBasis, SolverWorkspace, generate_basis, eval_bspline_kernel!
# Integrals
export init_scf_workspace, build_total_J_matrix!, assemble_K_matrix!, assemble_J_matrix!, build_specific_J_matrix!, build_specific_K_matrix!
# Poisson
export solve_poisson_J!, solve_generalized_poisson!
# Helpers
export evaluate_orbital
# CI
export extract_virtuals, get_cached_Rk!,compute_Rk

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

# CI
include("ci.jl")

end
