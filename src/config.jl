using LinearAlgebra

mutable struct Shell
    n::Int
    l::Int
    occ::Float64
    coeffs::Vector{Float64}
end

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
