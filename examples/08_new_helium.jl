using Pkg
Pkg.activate(joinpath(@__DIR__, "..")) 

using AtomicSplines
using LinearAlgebra
using Printf

function run_helium_scf()
    Z = 2.0
    helium = Atom(Z, [Orbital(1, 0, 2.0)])
    orb1s = helium.orbitals[1]

    # Basis Setup
    R_MAX = 20.0  # Increased slightly for better tail resolution
    N_ELEMS = 100
    ORDER = 6     # k=6 (Order 6 splines)
    basis = generate_basis(R_MAX, N_ELEMS, ORDER, γ=3.0) # Slightly sharper gamma

    # 1. Static Matrices
    T, V_nuc, S = assemble_core(basis, Z)
    H_core = T + V_nuc
    
    # Pre-calculate 1/r^2 for Poisson
    V_inv_r2 = assemble_operator_matrix(basis, r -> 1.0/r^2)

    # ======================================================================
    # CRITICAL FIX: Define Inner Indices
    # We solve only for coefficients 2 to n-1. 
    # c_1=0 fixes P(0)=0. c_n=0 fixes P(R)=0.
    # ======================================================================
    n = basis.num_splines
    inner = 2:(n-1)

    # 2. Initial Guess (He+)
    # Solve H c = E S c ONLY on the inner block
    H_dense = Matrix(H_core[inner, inner])
    S_dense = Matrix(S[inner, inner])

    vals, vecs = eigen(Symmetric(H_dense), Symmetric(S_dense))

   
    # Reconstruct full vector (padding zeros at boundaries)
    c_full = zeros(n)
    c_full[inner] = vecs[:, 1]
    
    orb1s.coeffs = c_full
    orb1s.energy = vals[1]

    println("Initial Guess (He+): E = $(orb1s.energy) (Should be ~ -2.0)")

    # 3. SCF Loop
    MAX_ITER = 40
    TOL = 1e-9
    MIXING = 0.6 # Increased mixing for stability
    E_old = 0.0

    println("\nIter | Orbital Energy (Ha) | Total Energy (Ha) | Delta")
    println("-"^60)

    for iter in 1:MAX_ITER
        old_coeffs = copy(orb1s.coeffs)

        # A. Poisson Solve (Returns full vector, BCs already handled inside)
        y_coeffs = solve_poisson_general(basis, orb1s.coeffs, orb1s.coeffs, 0, T, V_inv_r2)

        # B. Interaction Matrix
        V_ee = assemble_interaction_matrix(basis, y_coeffs)

        # C. Fock Matrix
        F = H_core + V_ee

        # D. Solve Eigenproblem (AGAIN: Only on inner block)

        F_dense = Matrix(F[inner, inner])
        # S_dense was computed above, reuse it

        vals, vecs = eigen(Symmetric(F_dense), Symmetric(S_dense))

        
        # Map back to full vector
        new_coeffs_inner = vecs[:, 1]
        
        # E. Mixing (Mix inner coeffs only, or full vector since boundaries are 0)
        c_full = zeros(n)
        c_full[inner] = new_coeffs_inner
        
        orb1s.energy = vals[1]
        orb1s.coeffs = MIXING .* c_full + (1 - MIXING) .* old_coeffs
        
        # Normalize (Always using full S matrix for safety)
        norm_val = sqrt(dot(orb1s.coeffs, S * orb1s.coeffs))
        orb1s.coeffs ./= norm_val

        # F. Total Energy
        V_ee_val = dot(orb1s.coeffs, V_ee * orb1s.coeffs)
        E_total = 2 * orb1s.energy - V_ee_val

        diff = abs(E_total - E_old)
        @printf("%4d | %18.10f | %17.10f | %8.2e\n", iter, orb1s.energy, E_total, diff)

        if diff < TOL
            println("-"^60)
            println("SCF Converged!")
            println("Final Total Energy: $(E_total) Ha")
            println("Reference He Limit: -2.86168 Ha")
            break
        end
        E_old = E_total
    end
end

run_helium_scf()
