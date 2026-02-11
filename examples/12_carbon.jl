using Pkg
Pkg.activate(joinpath(@__DIR__, "..")) 


using AtomicSplines
using LinearAlgebra
using Printf

# --- 1. The Projection Operator Helper ---
function assemble_projection_operator(basis, core_orbitals::Vector{Orbital}, S_matrix; eta=1e8)
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

function run_carbon_average_hf()
    println("==========================================")
    println("   CARBON (Z=6) - AVERAGE ENERGY HF")
    println("   Config: [He] 2s2 2p2")
    println("==========================================")

    # 1. Setup
    Z = 6.0; R_MAX = 8.0 # Carbon is smaller than Be
    # Higher order k=9 often helps with p-orbitals, but 7 is fine
    basis = generate_basis(R_MAX, 150, 7, γ=3.0) 
    
    n = basis.num_splines
    inner = 2:n  # Correct BCs for potentials (Neumann at R_max)
    inner_orb = 2:(n-1) # Dirichlet for orbitals

    # Operators
    T, V_nuc, S = assemble_core(basis, Z)
    V_inv_r2 = assemble_operator_matrix(basis, r -> 1.0/r^2)
    
    # H_core for s-orbitals (l=0)
    H_s = T + V_nuc
    
    # H_core for p-orbitals (l=1) -> Add 1/r^2 term
    # Centrifugal: l(l+1)/2 * 1/r^2. For l=1, factor is 1.0.
    H_p = H_s + 1.0 .* V_inv_r2

    # 2. Initial Guess (Hydrogenic / Screened)
    H_dense_s = Matrix(H_s[inner_orb, inner_orb])
    H_dense_p = Matrix(H_p[inner_orb, inner_orb])
    S_dense   = Matrix(S[inner_orb, inner_orb])

    # Diagonalize to get guess orbitals
    vals_s, vecs_s = eigen(Symmetric(H_dense_s), Symmetric(S_dense))
    vals_p, vecs_p = eigen(Symmetric(H_dense_p), Symmetric(S_dense))

    # Construct Initial Orbitals
    # 1s (Ground State S)
    c_1s = zeros(n); c_1s[inner_orb] = vecs_s[:, 1]
    orb1s = Orbital(1, 0, 2.0, vals_s[1], c_1s)

    # 2s (First Excited S)
    c_2s = zeros(n); c_2s[inner_orb] = vecs_s[:, 2]
    orb2s = Orbital(2, 0, 2.0, vals_s[2], c_2s)

    # 2p (Ground State P)
    c_2p = zeros(n); c_2p[inner_orb] = vecs_p[:, 1]
    orb2p = Orbital(2, 1, 2.0, vals_p[1], c_2p) # occ=2 for Carbon

    # 3. SCF Loop
    MAX_ITER = 40; TOL = 1e-7; MIXING = 0.6
    E_old = 0.0

    println("\nIter |  E(1s)  |  E(2s)  |  E(2p)  | Total Energy (Ha) | Delta")
    println("-"^80)

    for iter in 1:MAX_ITER
        # Save old for mixing
        c1s_old = copy(orb1s.coeffs)
        c2s_old = copy(orb2s.coeffs)
        c2p_old = copy(orb2p.coeffs)

        # --- A. POTENTIALS (Poisson) ---
        
        # 1. Monopoles (k=0) - Charge density
        y0_1s = solve_poisson_general(basis, orb1s.coeffs, orb1s.coeffs, 0, T, V_inv_r2)
        y0_2s = solve_poisson_general(basis, orb2s.coeffs, orb2s.coeffs, 0, T, V_inv_r2)
        y0_2p = solve_poisson_general(basis, orb2p.coeffs, orb2p.coeffs, 0, T, V_inv_r2)

        # 2. Exchange Dipoles (k=1) - For s-p interactions
        y1_1s_2p = solve_poisson_general(basis, orb1s.coeffs, orb2p.coeffs, 1, T, V_inv_r2)
        y1_2s_2p = solve_poisson_general(basis, orb2s.coeffs, orb2p.coeffs, 1, T, V_inv_r2)

        # 3. Self-Interaction Quadrupole (k=2) - For p-p interaction
        y2_2p_2p = solve_poisson_general(basis, orb2p.coeffs, orb2p.coeffs, 2, T, V_inv_r2)

        # --- B. MATRICES ---
        
        # Direct (J) - All spherical (k=0)
        J_1s = assemble_interaction_matrix(basis, y0_1s)
        J_2s = assemble_interaction_matrix(basis, y0_2s)
        J_2p = assemble_interaction_matrix(basis, y0_2p) # Average spherical part

        # Exchange (K)
        # s-s exchange (k=0)
        K_1s_s = assemble_exchange_matrix(basis, orb1s.coeffs, 0, T, V_inv_r2)
        K_2s_s = assemble_exchange_matrix(basis, orb2s.coeffs, 0, T, V_inv_r2)
        
        # s-p exchange (k=1)
        K_1s_p = assemble_exchange_matrix(basis, orb1s.coeffs, 1, T, V_inv_r2)
        K_2s_p = assemble_exchange_matrix(basis, orb2s.coeffs, 1, T, V_inv_r2)
        
        # p-s exchange (k=1) -> Same operator form, different orbital
        K_2p_to_s = assemble_exchange_matrix(basis, orb2p.coeffs, 1, T, V_inv_r2)

        # p-p self-exchange/shape (k=2)
        K_2p_shape = assemble_interaction_matrix(basis, y2_2p_2p) # Treated as potential for self

        # --- C. FOCK CONSTRUCTION (Average Energy Coefficients) ---
        
        # Shells: 1s(2), 2s(2), 2p(2)
        
        # F_1s (s-type)
        # H + J_1s + (2J_2s - K_2s) + (2J_2p - K_2p_avg)
        # Coeff for K(s,p) in avg energy is -1/3 * q_p
        # q_p = 2. So term is -2/3 K_2p_to_s
        F_1s = H_s + J_1s + (2.0 * J_2s - K_2s_s) + (2.0 * J_2p - (1.0/3.0) * K_2p_to_s)
        
        # F_2s (s-type)
        # H + (2J_1s - K_1s) + J_2s + (2J_2p - K_2p_avg)
        F_2s = H_s + (2.0 * J_1s - K_1s_s) + J_2s + (2.0 * J_2p - (1.0/3.0) * K_2p_to_s)
        
        # F_2p (p-type)
        # H_p + (2J_1s - K_1s_avg) + (2J_2s - K_2s_avg) + Self_Interaction
        # Coeff for K(p,s) is -1/3 * q_s = -2/3
        # Self Interaction for p^2 avg: 1.0 * J_2p - (2/25) * F2 term
        F_2p = H_p + (2.0 * J_1s - (1.0/3.0) * K_1s_p) + 
                     (2.0 * J_2s - (1.0/3.0) * K_2s_p) + 
                     (1.0 * J_2p - (2.0/25.0) * K_2p_shape)

        # --- D. SOLVE ---
        
        # 1. Solve 1s
        vals, vecs = eigen(Symmetric(Matrix(F_1s[inner_orb,inner_orb])), Symmetric(S_dense))
        c_1s_new = zeros(n); c_1s_new[inner_orb] = vecs[:, 1]; e_1s = vals[1]

        # 2. Solve 2s (Project against 1s)
        orb1s_temp = Orbital(1, 0, 2.0, e_1s, c_1s_new)
        P_op = assemble_projection_operator(basis, [orb1s_temp], S; eta=1e8)
        F_2s_proj = F_2s + P_op
        
        vals, vecs = eigen(Symmetric(Matrix(F_2s_proj[inner_orb,inner_orb])), Symmetric(S_dense))
        c_2s_new = zeros(n); c_2s_new[inner_orb] = vecs[:, 1]; e_2s = vals[1]

        # 3. Solve 2p (NO Projection needed! Orthogonal by symmetry)
        vals, vecs = eigen(Symmetric(Matrix(F_2p[inner_orb,inner_orb])), Symmetric(S_dense))
        c_2p_new = zeros(n); c_2p_new[inner_orb] = vecs[:, 1]; e_2p = vals[1]

        # --- E. UPDATE & ORTHOGONALIZE ---
        
        # Mix
        orb1s.coeffs = MIXING * c_1s_new + (1-MIXING) * c1s_old
        orb2s.coeffs = MIXING * c_2s_new + (1-MIXING) * c2s_old
        orb2p.coeffs = MIXING * c_2p_new + (1-MIXING) * c2p_old
        
        # Normalize & Orthogonalize (Gram-Schmidt for s-orbitals)
        orb1s.coeffs ./= sqrt(dot(orb1s.coeffs, S * orb1s.coeffs))
        
        # 2s against 1s
        overlap = dot(orb1s.coeffs, S * orb2s.coeffs)
        orb2s.coeffs .-= overlap .* orb1s.coeffs
        orb2s.coeffs ./= sqrt(dot(orb2s.coeffs, S * orb2s.coeffs))
        
        # 2p (Just Normalize)
        orb2p.coeffs ./= sqrt(dot(orb2p.coeffs, S * orb2p.coeffs))
        
        # Store Energies
        orb1s.energy = e_1s; orb2s.energy = e_2s; orb2p.energy = e_2p

        # --- F. TOTAL ENERGY ---
        # E_avg = Sum(q*h) + Sum(Pairs)
        # Note: This is simpler than summing F matrix elements which double count
        
        # One-electron terms
        h1s = dot(orb1s.coeffs, H_s * orb1s.coeffs)
        h2s = dot(orb2s.coeffs, H_s * orb2s.coeffs)
        h2p = dot(orb2p.coeffs, H_p * orb2p.coeffs)
        E_one = 2*h1s + 2*h2s + 2*h2p
        
        # Two-electron terms (Average Energy Coeffs)
        # 1s-1s: F0
        E_1s1s = 1.0 * dot(orb1s.coeffs, J_1s * orb1s.coeffs)
        # 2s-2s: F0
        E_2s2s = 1.0 * dot(orb2s.coeffs, J_2s * orb2s.coeffs)
        # 2p-2p: F0 - (2/25)F2
        F0_pp = dot(orb2p.coeffs, J_2p * orb2p.coeffs)
        F2_pp = dot(orb2p.coeffs, K_2p_shape * orb2p.coeffs)
        E_2p2p = 1.0 * F0_pp - (2.0/25.0) * F2_pp
        
        # Inter-shell
        # 1s-2s: 4*F0 - 2*G0
        J_1s2s = dot(orb1s.coeffs, J_2s * orb1s.coeffs)
        K_1s2s = dot(orb1s.coeffs, K_2s_s * orb1s.coeffs)
        E_1s2s = 4*J_1s2s - 2*K_1s2s
        
        # 1s-2p: 4*F0 - (2/3)*G1
        J_1s2p = dot(orb1s.coeffs, J_2p * orb1s.coeffs)
        # Careful: G1 = <1s|K_2p(k=1)|1s>
        G1_1s2p = dot(orb1s.coeffs, K_2p_to_s * orb1s.coeffs)
        E_1s2p = 4*J_1s2p - (2.0/3.0)*G1_1s2p
        
        # 2s-2p: 4*F0 - (2/3)*G1
        J_2s2p = dot(orb2s.coeffs, J_2p * orb2s.coeffs)
        G1_2s2p = dot(orb2s.coeffs, K_2p_to_s * orb2s.coeffs)
        E_2s2p = 4*J_2s2p - (2.0/3.0)*G1_2s2p
        
        E_total = E_one + E_1s1s + E_2s2s + E_2p2p + E_1s2s + E_1s2p + E_2s2p
        
        diff = abs(E_total - E_old)
        @printf("%4d | %9.5f | %9.5f | %9.5f | %14.8f | %8.2e\n", 
                iter, e_1s, e_2s, e_2p, E_total, diff)
        
        if diff < TOL
            println("-"^80)
            println("Converged!")
            println("Ref (Free C Avg): ~ -37.688 Ha")
            break
        end
        E_old = E_total
    end
end

run_carbon_average_hf()
