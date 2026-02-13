using LinearAlgebra

# --- 1. Data Structures ---

mutable struct Shell
    n::Int
    l::Int
    occ::Float64
    coeffs::Vector{Float64} # Default constructor expects (Int, Int, Float64, Vector)
end

# FIXED CONSTRUCTOR:
# We explicitly specify 'basis_size::Int' so this signature 
# (Int, Int, Float64, Int) is distinct from the default (..., Vector).
Shell(n::Int, l::Int, occ::Float64, basis_size::Int) = Shell(n, l, occ, zeros(basis_size))

struct InteractionTerm
    k::Int
    coeff::Float64
    is_exchange::Bool
end

struct AtomicConfig
    Z::Float64
    shells::Vector{Shell}
end
