export Orbital, Atom
export solve_eigen

# The Orbital Struct
mutable struct Orbital
    n::Int          # Principal quantum number
    l::Int          # Angular momentum
    occ::Float64    # Occupation
    energy::Float64 # Eigenenergy
    coeffs::Vector{Float64} # B-spline coefficients
end

# Constructor for an empty/initial orbital
Orbital(n, l, occ) = Orbital(n, l, occ, 0.0, Float64[])

# The Atom Struct
struct Atom
    Z::Float64
    orbitals::Vector{Orbital}
end

Atom(Z) = Atom(Z, Orbital[])

# ===================================================================
# Quantum State Data Structures
# ===================================================================

"""
Immutable structure representing an atomic LS Term.
L: Total orbital angular momentum quantum number (0=S, 1=P, 2=D, etc.)
two_S: Twice the total spin quantum number (2S) to maintain integer exactness.
"""
struct LSTerm
    L::Int
    two_S::Int
end
