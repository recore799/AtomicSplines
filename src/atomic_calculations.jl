# ===================================================================
# Fractional Parentage and Tensor Scaling Logic
# ===================================================================

# --- Pre-defined Basis States ---

# Parent Configuration States (p^2)
const parent_3P = LSTerm(1, 2)
const parent_1D = LSTerm(2, 0)
const parent_1S = LSTerm(0, 0)

# Child Configuration States (p^3)
const child_4S  = LSTerm(0, 3)
const child_2D  = LSTerm(2, 1)
const child_2P  = LSTerm(1, 1)

# --- Theoretical Tabulations ---

"""
Associative mapping from a p^3 child term to an iterable collection of 
its allowed p^2 parent terms and their exact CFP amplitudes.
"""
const CFP_p3 = Dict{LSTerm, Vector{Tuple{LSTerm, Float64}}}(
    child_4S => [
        (parent_3P, 1.0)
    ],
    child_2D => [
        (parent_3P, -1.0 / sqrt(2.0)),
        (parent_1D,  1.0 / sqrt(2.0))
    ],
    child_2P => [
        (parent_3P, -1.0 / sqrt(2.0)),
        (parent_1D,  sqrt(5.0 / 18.0)),
        (parent_1S,  sqrt(2.0 / 9.0))
    ]
)

"""
Associative mapping from a p^2 parent term to its exact theoretical 
energy shift relative to the configuration average (E_av) for F^2.
"""
const parent_F2_shifts = Dict{LSTerm, Rational{Int}}(
    parent_3P => -5 // 25,
    parent_1D =>  1 // 25,
    parent_1S => 10 // 25
)

# --- Dynamic Scaling Algorithms ---

"""
    get_parent_amplitudes(target_term::LSTerm)

Retrieves the fractional parentage pathways for a specified term.
"""
function get_parent_amplitudes(target_term::LSTerm)
    if haskey(CFP_p3, target_term)
        return CFP_p3[target_term]
    else
        error("Quantum State Error: Unmapped term requested.")
    end
end

"""
    compute_child_F2_coefficient(child_term::LSTerm, n::Int)

Computes the absolute F^2 exchange coefficient for an l^n configuration
by summing the exact fractional contributions of its parent states.
Returns the coefficient as a Float64 for immediate use in matrix assembly.
"""
function compute_child_F2_coefficient(child_term::LSTerm, n::Int)
    if n <= 2
        error("Scaling Theorem Error: The prefactor n/(n-2) is invalid for n <= 2.")
    end
    
    prefactor = n // (n - 2) 
    parent_pathways = get_parent_amplitudes(child_term)
    
    total_F2_coefficient = 0 // 1
    
    for (parent, amplitude) in parent_pathways
        parent_shift = parent_F2_shifts[parent]
        
        # Isolate the exact rational amplitude squared to prevent truncation
        sq_amp_rational = rationalize(amplitude^2, tol=1e-8) 
        total_F2_coefficient += sq_amp_rational * parent_shift
    end
    
    exact_fraction = prefactor * total_F2_coefficient
    
    # Cast to Float64 for direct application to the B-spline integration matrices
    return Float64(exact_fraction)
end

function compute_effective_central_potential(ws::SolverWorkspace, orbitals::Vector{Orbital}, dense_grid::Vector{Float64}, Z::Float64)
    n_splines = ws.basis.num_splines
    total_Y0_coeffs = zeros(Float64, n_splines)
    y_temp = zeros(Float64, n_splines)
    
    # 1. Superposition of the Hartree Screening Functions
    for orb in orbitals
        if orb.occ > 0.0
            # Solve the Poisson equation for the individual orbital's density
            solve_poisson_J!(ws, y_temp, orb)
            
            # Linearly scale by fractional occupancy and accumulate
            total_Y0_coeffs .+= orb.occ .* y_temp
        end
    end
    
    # 2. Evaluate the continuous spline over the dense physical grid
    Y0_grid = evaluate_orbital(ws.basis, total_Y0_coeffs, dense_grid)
    
    # 3. Construct V_eff(r) = (-Z + Y0(r)) / r
    n_points = length(dense_grid)
    V_eff = zeros(Float64, n_points)
    
    for i in 1:n_points
        r = dense_grid[i]
        # Prevent division by zero precisely at the origin
        if r > 1e-12
            V_eff[i] = (-Z + Y0_grid[i]) / r
        else
            # Deep in the core, the potential is heavily dominated by the bare nucleus
            V_eff[i] = -Z / 1e-12 
        end
    end
    
    return V_eff
end

