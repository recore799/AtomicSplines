module AtomicSplines

using LinearAlgebra
using BandedMatrices
using FastGaussQuadrature
using StaticArrays
using JLD2
using Printf

# Basis
export BSplineBasis, SolverWorkspace, generate_basis, eval_bspline_kernel!
# Integrals
export init_scf_workspace, cached_init_scf_workspace, build_total_J_matrix!, assemble_K_matrix!, assemble_J_matrix!, build_specific_J_matrix!, build_specific_K_matrix!
# Poisson
export solve_poisson_J!, solve_generalized_poisson!
# Helpers
export evaluate_orbital
# --- CFP ---
export LSTerm, get_parent_amplitudes, compute_child_F2_coefficient

# CI
export extract_virtuals, get_cached_Rk!,compute_Rk

# SO
export compute_effective_central_potential

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

# High level logic
include("atomic_calculations.jl")
include("ci.jl")

end
