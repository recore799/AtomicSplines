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

# --- 2. Main Beryllium Routine ---
function run_beryllium_full_hf()
    println("==========================================")
    println("   BERYLLIUM (Z=4) - FULL HARTREE-FOCK")
    println("   (Direct + Exchange + Projection)")
    println("==========================================")

    # 1. Setup (Same as before)
    Z = 4.0; R_MAX = 2.0
    basis = generate_basis(R_MAX, 150, 7, γ=2.5)
    T, V_nuc, S = assemble_core(basis, Z)
    H_core = T + V_nuc
    V_inv_r2 = assemble_operator_matrix(basis, r -> 1.0/r^2)
    n = basis.num_splines; inner = 2:(n-1)

    # 2. Initial Guess (From previous Hartree result or He-like)
    # Let's do a quick He-like guess to keep it self-contained
    H_dense = Matrix(H_core[inner, inner]); S_dense = Matrix(S[inner, inner])
    vals, vecs = eigen(Symmetric(H_dense), Symmetric(S_dense))
    
    c_1s = zeros(n); c_1s[inner] = vecs[:, 1]
    orb1s = Orbital(1, 0, 2.0, vals[1], c_1s)
    
    c_2s = zeros(n); c_2s[inner] = vecs[:, 2] 
    orb2s = Orbital(2, 0, 2.0, vals[2], c_2s)

    # 3. SCF Loop
    MAX_ITER = 30; TOL = 1e-8; MIXING = 1.0; E_old = 0.0
    
    println("\nIter |   E(1s)   |   E(2s)   | Total Energy (Ha) | Delta")
    println("-"^70)

    for iter in 1:MAX_ITER
        # Save old for mixing
        c1s_old = copy(orb1s.coeffs); c2s_old = copy(orb2s.coeffs)

        # --- A. Build Direct Matrices (J) ---
        y_1s = solve_poisson_general(basis, orb1s.coeffs, orb1s.coeffs, 0, T, V_inv_r2)
        J_1s = assemble_interaction_matrix(basis, y_1s)

        y_2s = solve_poisson_general(basis, orb2s.coeffs, orb2s.coeffs, 0, T, V_inv_r2)
        J_2s = assemble_interaction_matrix(basis, y_2s)

        # --- B. Build Exchange Matrices (K) ---
        # Exchange between 1s and 2s (k=0 because both are s-orbitals)
        # K_1s_op: The operator K_{1s} acting on other functions
        K_1s_op = assemble_exchange_matrix(basis, orb1s.coeffs, 0, T, V_inv_r2)
        
        # K_2s_op: The operator K_{2s} acting on other functions
        K_2s_op = assemble_exchange_matrix(basis, orb2s.coeffs, 0, T, V_inv_r2)

        # --- C. Solve 1s Orbital ---
        # F_1s = H + J_1s + (2J_2s - K_2s)
        F_1s = H_core + J_1s + (2.0 * J_2s - K_2s_op)
        
        vals_1s, vecs_1s = eigen(Symmetric(Matrix(F_1s[inner,inner])), Symmetric(S_dense))
        c_1s_new = zeros(n); c_1s_new[inner] = vecs_1s[:, 1]
        e_1s_new = vals_1s[1]

        # --- D. Solve 2s Orbital (With Projector) ---
        # F_2s = H + (2J_1s - K_1s) + J_2s
        F_2s = H_core + (2.0 * J_1s - K_1s_op) + J_2s
        
        # Project against the NEW 1s (better stability)
        orb1s_temp = Orbital(1, 0, 2.0, e_1s_new, c_1s_new)
        P_op = assemble_projection_operator(basis, [orb1s_temp], S; eta=1000.0)
        
        F_2s_proj = F_2s + P_op
        
        vals_2s, vecs_2s = eigen(Symmetric(Matrix(F_2s_proj[inner,inner])), Symmetric(S_dense))
        c_2s_new = zeros(n); c_2s_new[inner] = vecs_2s[:, 1]
        e_2s_new = vals_2s[1]

        # --- E. Mix & Normalize ---
        orb1s.coeffs = MIXING * c_1s_new + (1-MIXING) * c1s_old
        orb2s.coeffs = MIXING * c_2s_new + (1-MIXING) * c2s_old
        
        orb1s.coeffs ./= sqrt(dot(orb1s.coeffs, S * orb1s.coeffs))
        orb2s.coeffs ./= sqrt(dot(orb2s.coeffs, S * orb2s.coeffs))
        
        orb1s.energy = e_1s_new
        orb2s.energy = e_2s_new

        # --- F. Total Energy (Full HF) ---
        # E = 2*h1 + 2*h2 + F11 + F22 + 4*F12 - 2*G12
        
        h1 = dot(orb1s.coeffs, H_core * orb1s.coeffs)
        h2 = dot(orb2s.coeffs, H_core * orb2s.coeffs)
        
        f11 = dot(orb1s.coeffs, J_1s * orb1s.coeffs)
        f22 = dot(orb2s.coeffs, J_2s * orb2s.coeffs)
        f12 = dot(orb1s.coeffs, J_2s * orb1s.coeffs)
        
        # Exchange Scalar G12 = <1s | K_2s | 1s>
        g12 = dot(orb1s.coeffs, K_2s_op * orb1s.coeffs)
        
        E_total = 2*h1 + 2*h2 + f11 + f22 + 4*f12 - 2*g12
        
        diff = abs(E_total - E_old)
        @printf("%4d | %9.5f | %9.5f | %17.10f | %8.2e\n", 
                iter, orb1s.energy, orb2s.energy, E_total, diff)
        
        if diff < TOL
            println("-"^70)
            println("Converged!")
            println("Final Energy: $E_total Ha")
            println("Target Limit: -14.573023 Ha")
            break
        end
        E_old = E_total
    end
end

run_beryllium_full_hf()


