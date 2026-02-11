using Pkg
Pkg.activate(joinpath(@__DIR__, "..")) 


using AtomicSplines
using LinearAlgebra
using Printf

# --- 1. The Projection Operator Helper ---
function assemble_projection_operator(basis, core_orbitals::Vector{Orbital}, S_matrix; eta=1e6)
    n = basis.num_splines
    P_mat = zeros(n, n)
    
    for orb in core_orbitals
        # P = eta * |core><core|
        # In generalized basis: P = eta * (S * c) * (c' * S)
        ket = S_matrix * orb.coeffs
        P_mat .+= eta .* (ket * ket')
    end
    
    return P_mat
end

# --- 2. Main Beryllium Routine ---
function run_beryllium_hartree()
    println("==========================================")
    println("   BERYLLIUM (Z=4) - HARTREE APPROX")
    println("   (Orthogonality via Projection)")
    println("==========================================")

    # A. Setup Basis
    Z = 4.0
    R_MAX = 20.0
    basis = generate_basis(R_MAX, 120, 7, γ=2.5) # Order 7 for high accuracy
    
    # B. Static Matrices
    T, V_nuc, S = assemble_core(basis, Z)
    H_core = T + V_nuc
    V_inv_r2 = assemble_operator_matrix(basis, r -> 1.0/r^2)
    
    # Slice helper for Boundary Conditions
    n = basis.num_splines
    inner = 2:(n-1)

    # C. Initial Guess
    # 1. Solve He-like 1s (Z=4) ignoring 2s
    #    (Standard dense solve on inner block)
    H_dense = Matrix(H_core[inner, inner])
    S_dense = Matrix(S[inner, inner])
    vals, vecs = eigen(Symmetric(H_dense), Symmetric(S_dense))
    
    c_1s = zeros(n)
    c_1s[inner] = vecs[:, 1]
    orb1s = Orbital(1, 0, 2.0, vals[1], c_1s)
    
    # 2. Guess 2s
    #    We can just take the 2nd eigenvector of the Z=4 system as a rough start
    #    It's physically wrong (doesn't see screening), but it's orthogonal to 1s.
    c_2s = zeros(n)
    c_2s[inner] = vecs[:, 2] 
    orb2s = Orbital(2, 0, 2.0, vals[2], c_2s)
    
    println("Initial Guess Prepared.")

    # D. SCF Loop
    MAX_ITER = 50
    TOL = 1e-8
    MIXING = 0.5
    E_old = 0.0
    
    println("\nIter |   E(1s)   |   E(2s)   | Total Energy (Ha) | Delta")
    println("-"^70)

    for iter in 1:MAX_ITER
        # Store old coefficients for mixing
        c1s_old = copy(orb1s.coeffs)
        c2s_old = copy(orb2s.coeffs)

        # ------------------------------------------------------------------
        # STEP 1: Build Potentials
        # ------------------------------------------------------------------
        
        # Potential Y(1s, 1s) -> Seen by 2s (strength 2) and 1s (strength 1)
        y_1s = solve_poisson_general(basis, orb1s.coeffs, orb1s.coeffs, 0, T, V_inv_r2)
        J_1s_mat = assemble_interaction_matrix(basis, y_1s)

        # Potential Y(2s, 2s) -> Seen by 1s (strength 2) and 2s (strength 1)
        y_2s = solve_poisson_general(basis, orb2s.coeffs, orb2s.coeffs, 0, T, V_inv_r2)
        J_2s_mat = assemble_interaction_matrix(basis, y_2s)

        # ------------------------------------------------------------------
        # STEP 2: Solve 1s Orbital
        # ------------------------------------------------------------------
        # 1s feels: Core + 1 other 1s + 2 2s electrons
        F_1s = H_core + J_1s_mat + 2.0 * J_2s_mat
        
        # Solve Dense
        F_1s_dense = Matrix(F_1s[inner, inner])
        vals_1s, vecs_1s = eigen(Symmetric(F_1s_dense), Symmetric(S_dense))
        
        # Update 1s
        c_1s_new = zeros(n)
        c_1s_new[inner] = vecs_1s[:, 1]
        e_1s_new = vals_1s[1]

        # ------------------------------------------------------------------
        # STEP 3: Solve 2s Orbital (WITH PROJECTION)
        # ------------------------------------------------------------------
        # 2s feels: Core + 2 1s electrons + 1 other 2s electron
        F_2s = H_core + 2.0 * J_1s_mat + J_2s_mat
        
        # Build Projector: P = |1s_new><1s_new| (Using the updated 1s is better)
        # Note: Need updated orbital struct for helper
        orb1s_temp = Orbital(1, 0, 2.0, e_1s_new, c_1s_new)
        P_op = assemble_projection_operator(basis, [orb1s_temp], S; eta=100.0) # eta=100 is usually enough
        
        F_2s_proj = F_2s + P_op
        
        # Solve Dense
        F_2s_dense = Matrix(F_2s_proj[inner, inner])
        vals_2s, vecs_2s = eigen(Symmetric(F_2s_dense), Symmetric(S_dense))
        
        # Update 2s (Take lowest root! Projector moved 1s to ~100)
        c_2s_new = zeros(n)
        c_2s_new[inner] = vecs_2s[:, 1]
        e_2s_new = vals_2s[1]

        # ------------------------------------------------------------------
        # STEP 4: Update & Mix
        # ------------------------------------------------------------------
        # Mix
        orb1s.coeffs = MIXING * c_1s_new + (1-MIXING) * c1s_old
        orb2s.coeffs = MIXING * c_2s_new + (1-MIXING) * c2s_old
        
        # Normalize
        orb1s.coeffs ./= sqrt(dot(orb1s.coeffs, S * orb1s.coeffs))
        orb2s.coeffs ./= sqrt(dot(orb2s.coeffs, S * orb2s.coeffs))
        
        orb1s.energy = e_1s_new
        orb2s.energy = e_2s_new

        # ------------------------------------------------------------------
        # STEP 5: Calculate Total Energy (Hartree)
        # ------------------------------------------------------------------
        # E = 2*h_1s + 2*h_2s + F(1s,1s) + F(2s,2s) + 4*F(1s,2s)
        
        # One-electron terms <i|h|i>
        h1 = dot(orb1s.coeffs, H_core * orb1s.coeffs)
        h2 = dot(orb2s.coeffs, H_core * orb2s.coeffs)
        
        # Direct Integrals F(a,b) = <a|J_b|a>
        # Note: We use the matrices J calculated from current density
        # Since we mixed coefficients, strictly we should recompute J, 
        # but using the ones from this step is fine for convergence check.
        f11 = dot(orb1s.coeffs, J_1s_mat * orb1s.coeffs)
        f22 = dot(orb2s.coeffs, J_2s_mat * orb2s.coeffs)
        f12 = dot(orb1s.coeffs, J_2s_mat * orb1s.coeffs) # Symmetric
        
        E_total = 2*h1 + 2*h2 + f11 + f22 + 4*f12
        
        diff = abs(E_total - E_old)
        @printf("%4d | %9.5f | %9.5f | %17.10f | %8.2e\n", 
                iter, orb1s.energy, orb2s.energy, E_total, diff)
        
        if diff < TOL
            println("-"^70)
            println("Converged!")
            println("Total Energy: $E_total Ha")
            println("Virial Ratio: ... (todo)")
            break
        end
        E_old = E_total
    end
end

run_beryllium_hartree()
