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
