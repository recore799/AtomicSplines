using Pkg
Pkg.activate(joinpath(@__DIR__, "..")) 

using AtomicSplines
using LinearAlgebra
using SparseArrays
using Printf
using Dates # Added for timing

 function run_neon()

    println("==========================================")

    println("      NEON (Z=10) - AUTOMATED HF")

    println("==========================================")


    # 1. Setup

    Z = 10.0

    R_MAX = 10.0

    basis = generate_basis(R_MAX, 120, 7, γ=3.0) # More splines for p-orbitals

    n = basis.num_splines

    inner = 2:(n-1)


    # 2. Define Configuration

    # Neon: 1s2, 2s2, 2p6

    shells = [

        Shell(1, 0, 2.0, n),

        Shell(2, 0, 2.0, n),

        Shell(2, 1, 6.0, n)

    ]

    config = AtomicConfig(Z, shells)


    # 3. Static Operators

    T, V_nuc, S = assemble_core(basis, Z)

    V_inv_r2 = assemble_operator_matrix(basis, r -> 1.0/r^2)

    

    # Pre-compute H_core for l=0 and l=1

    H_core_l0 = T + V_nuc

    H_core_l1 = T + V_nuc + assemble_centrifugal(basis, 1)

    

    H_cores = Dict(0 => H_core_l0, 1 => H_core_l1)

    

    # 4. Initial Guess (Screened Hydrogenic is better than bare nucleus for Neon)

    # But for simplicity, let's just solve bare nucleus H_core

    # 1s & 2s (l=0)

    vals_s, vecs_s = eigen(Symmetric(Matrix(H_core_l0[inner,inner])), Symmetric(Matrix(S[inner,inner])))

    config.shells[1].coeffs[inner] = vecs_s[:, 1] # 1s

    config.shells[2].coeffs[inner] = vecs_s[:, 2] # 2s (Orthogonal to 1s by eigen)

    

    # 2p (l=1)

    vals_p, vecs_p = eigen(Symmetric(Matrix(H_core_l1[inner,inner])), Symmetric(Matrix(S[inner,inner])))

    config.shells[3].coeffs[inner] = vecs_p[:, 1] # 2p (Lowest l=1 state)

    

    println("Initial Guess Generated.")


    # 5. SCF Loop

    MAX_ITER = 40

    MIXING = 0.3

    E_old = 0.0


    println("\nIter | Time  |  E(1s)   |   E(2s)   |   E(2p)   | Total Energy | Delta")

    println("-"^75)


    for iter in 1:MAX_ITER

        t_start = time()

        # Store old coeffs

        old_coeffs = [copy(s.coeffs) for s in config.shells]


        # --- A. Pre-Compute Potentials ---

        # 1. Direct Matrices (J) for each shell

        J_matrices = []

        for s in config.shells

            # Solve Y^0(s, s)

            y_vec = solve_poisson_general(basis, s.coeffs, s.coeffs, 0, T, V_inv_r2)

            push!(J_matrices, assemble_interaction_matrix(basis, y_vec))

        end

        

        # 2. Exchange Accessor (Lazy evaluation or Pre-computed?)

        # Since we might need K^0, K^1, K^2... let's compute on demand but cache them per iter.

        # We'll use a Dictionary: (shell_idx, k) -> Matrix

        K_cache = Dict{Tuple{Int, Int}, SparseMatrixCSC{Float64, Int64}}()

        

        get_K = (source_idx, k) -> begin

            key = (source_idx, k)

            if !haskey(K_cache, key)

                coeffs = config.shells[source_idx].coeffs

                # Compute and store

                K_cache[key] = assemble_exchange_matrix(basis, coeffs, k, T, V_inv_r2)

            end

            return K_cache[key]

        end


        # --- B. Update Each Shell ---

        energies = zeros(3)

        new_coeff_list = []

        

        for (i, shell) in enumerate(config.shells)

            # 1. Build Fock Matrix

            F = assemble_fock_matrix(basis, config, i, H_cores, J_matrices, get_K)

            

            # 2. Apply Projection?

            # 2s (i=2) must be orthogonal to 1s (i=1)

            # 1s (i=1) needs no projection

            # 2p (i=3) needs no projection (l=1 is orthogonal to l=0)

            if i == 2

                # Project out 1s (using the UPDATED 1s if available? 

                # Standard HF uses current iter's best guess. Let's use old 1s for stability or new?)

                # Let's use old 1s to match standard "freeze" logic, or better:

                # If we updated 1s already, use new.

                # Since we are in loop, 1s is at index 1.

                # We haven't updated config.shells[1] yet in place. We have `new_coeff_list`.

                # Let's use the `old_coeffs[1]` for simplicity.

                

                # Note: pass a dummy orbital struct or vector

                # We can reuse our `assemble_projection_operator`.

                # We need to wrap the coeff in an Orbital-like object or modify the function.

                # Quick hack: Make a temporary orbital

                temp_orb = Orbital(1, 0, 0.0, 0.0, old_coeffs[1])

                P = assemble_projection_operator(basis, [temp_orb], S; eta=1000.0)

                F .+= P

            end

            

            # 3. Solve

            vals, vecs = eigen(Symmetric(Matrix(F[inner,inner])), Symmetric(Matrix(S[inner,inner])))

            

            # Store

            energies[i] = vals[1]

            c_new = zeros(n)

            c_new[inner] = vecs[:, 1]

            push!(new_coeff_list, c_new)

        end

        

        # --- C. Update State ---

        for i in 1:3

            # Mix

            config.shells[i].coeffs = MIXING * new_coeff_list[i] + (1-MIXING) * old_coeffs[i]

            # Normalize

            norm_val = sqrt(dot(config.shells[i].coeffs, S * config.shells[i].coeffs))

            config.shells[i].coeffs ./= norm_val

        end


        # --- D. Total Energy (Hardcoded for Neon) ---

        # Shells: 1=1s, 2=2s, 3=2p

        s1 = config.shells[1]; s2 = config.shells[2]; s3 = config.shells[3]

        

        # 1. One-Electron Energies (Kinetic + Nuclear)

        # H_cores[0] is for s-orbitals, H_cores[1] for p-orbitals

        h1 = dot(s1.coeffs, H_cores[0] * s1.coeffs)

        h2 = dot(s2.coeffs, H_cores[0] * s2.coeffs)

        h3 = dot(s3.coeffs, H_cores[1] * s3.coeffs)

        

        E_one = s1.occ * h1 + s2.occ * h2 + s3.occ * h3


        # 2. Interaction Energies (Direct F and Exchange G)

        # Notation: F(a,b) is interaction between full shells a and b

        

        # We need the scalars: <a | J_b | a> and <a | K_b | a>

        # We already computed J_matrices[i] which represents J_i (potential of shell i)

        

        # Direct Terms (F_ab)

        # Self-Direct: <1s|J1s|1s>, <2s|J2s|2s>, <2p|J2p|2p>

        F11 = dot(s1.coeffs, J_matrices[1] * s1.coeffs)

        F22 = dot(s2.coeffs, J_matrices[2] * s2.coeffs)

        F33 = dot(s3.coeffs, J_matrices[3] * s3.coeffs)

        

        # Cross-Direct: <1s|J2s|1s>, etc.

        F12 = dot(s1.coeffs, J_matrices[2] * s1.coeffs)

        F13 = dot(s1.coeffs, J_matrices[3] * s1.coeffs)

        F23 = dot(s2.coeffs, J_matrices[3] * s2.coeffs)


        # Exchange Terms (G_ab)

        # We need to explicitly compute the scalar <a|K_b|a>

        # NOTE: get_K(source, k) returns the matrix operator.

        

        # 1s-2s (k=0)

        K_1s_2s = dot(s1.coeffs, get_K(2, 0) * s1.coeffs)

        

        # 1s-2p (k=1) - Note: get_K(3, 1) gets K operator from 2p(source) with k=1

        K_1s_2p = dot(s1.coeffs, get_K(3, 1) * s1.coeffs)

        

        # 2s-2p (k=1)

        K_2s_2p = dot(s2.coeffs, get_K(3, 1) * s2.coeffs)

        

        # 2p-2p Self Exchange (k=0 and k=2)

        K_2p_self_0 = dot(s3.coeffs, get_K(3, 0) * s3.coeffs)

        K_2p_self_2 = dot(s3.coeffs, get_K(3, 2) * s3.coeffs)


        # 3. Summation (Coefficients from Slater integrals for closed shells)

        # E_int = Sum_{pairs} (Direct - Exchange)

        # Note: The factors below account for the number of electrons and permutations

        # For closed shells A and B:

        # E(A,A) = N_A * (N_A - 1) / 2 * F0 - ... (Complex!)

        

        # EASIER WAY:

        # Use the relation E_int = 0.5 * Sum_i N_i (V_i - K_i)

        # V_i = Pot seen by shell i

        

        # Let's use the orbital energies directly to simplify this messy sum:

        # E_total = Sum(N_i * eps_i) - E_ee

        # E_ee = 0.5 * Sum_i N_i * <i | V_HF_i - H_core | i >

        #      = 0.5 * Sum_i N_i * (eps_i - h_i)

        

        # This is strictly valid for closed shell HF!

        

        term1 = s1.occ * (energies[1] - h1)

        term2 = s2.occ * (energies[2] - h2)

        term3 = s3.occ * (energies[3] - h3)

        

        E_ee = 0.5 * (term1 + term2 + term3)

        

        E_total = E_one + E_ee

        # Update print
        t_end = time() # Stop timer
        elapsed = t_end - t_start
        
        diff = abs(E_total - E_old)
        
        @printf("%4d | %6.2fs | %9.5f | %9.5f | %9.5f | %12.8f | %8.2e\n", 
                iter, elapsed, energies[1], energies[2], energies[3], E_total, diff)

       


        if diff < 1e-7 && iter > 5

            println("Converged!")

            # Calculate Total Energy properly once at the end?

            break

        end

        E_old = energies[1]

    end

end 
run_neon()
