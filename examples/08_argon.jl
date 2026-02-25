using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using LinearAlgebra
using Printf

function solve_argon(R_max)
    println("=== Argon (Z=18) RHF Solver ===")
    
    # 1. Setup
    N_elems = 200  # Increased for the tighter core and more nodes
    Z = 18.0
    basis = generate_basis(R_max, N_elems, Val(7), γ=3.0)
    ws = init_scf_workspace(basis, Z)

    empty!(ws.poisson_factors)

    n = basis.num_splines

    # Active sets based on boundary condition P(r) ~ r^(l+1)
    active_s = 2:(n-1)  # Drops B1. Allows B2 (r^1)
    active_p = 3:(n-1)  # Drops B1, B2. Starts at B3 (r^2)
    
    # Initialize Orbitals (n, l, occ)
    orbitals = [
        Orbital(1, 0, 2.0), # 1s
        Orbital(2, 0, 2.0), # 2s
        Orbital(2, 1, 6.0), # 2p
        Orbital(3, 0, 2.0), # 3s
        Orbital(3, 1, 6.0)  # 3p
    ]
    
    # Core Hamiltonians
    H_core_s = ws.T + ws.V
    H_core_p = ws.T + ws.V + ws.V2 
    
    # Initial Guess
    evals_s, evecs_s = eigen(Symmetric(H_core_s[active_s, active_s]), Symmetric(ws.S[active_s, active_s]))
    orbitals[1].coeffs = zeros(Float64, n); orbitals[1].coeffs[active_s] = evecs_s[:, 1]
    orbitals[2].coeffs = zeros(Float64, n); orbitals[2].coeffs[active_s] = evecs_s[:, 2]
    orbitals[4].coeffs = zeros(Float64, n); orbitals[4].coeffs[active_s] = evecs_s[:, 3] # 3s
    
    evals_p, evecs_p = eigen(Symmetric(H_core_p[active_p, active_p]), Symmetric(ws.S[active_p, active_p]))
    orbitals[3].coeffs = zeros(Float64, n); orbitals[3].coeffs[active_p] = evecs_p[:, 1]
    orbitals[5].coeffs = zeros(Float64, n); orbitals[5].coeffs[active_p] = evecs_p[:, 2] # 3p
    
    for orb in orbitals
        orb.coeffs ./= sqrt(dot(orb.coeffs, ws.S * orb.coeffs))
    end

    E_old = 0.0
    MIXING = 0.0
    
    println("Starting SCF Loop...")
    
    for iter in 1:100
        t0 = time()
        
        # --- BUILD POTENTIALS ---
        build_total_J_matrix!(ws, orbitals)
        
        assemble_K_matrix!(ws, ws.K_mats[0], 0, orbitals)
        assemble_K_matrix!(ws, ws.K_mats[1], 1, orbitals)
        
        # --- BUILD FOCK MATRICES ---
        ws.F_s .= H_core_s .+ ws.J .- ws.K_mats[0]
        ws.F_p .= H_core_p .+ ws.J .- ws.K_mats[1]
        
        # --- DIAGONALIZE ---
        evals_fs, evecs_fs = eigen(Symmetric(ws.F_s[active_s, active_s]), Symmetric(ws.S[active_s, active_s]))
        evals_fp, evecs_fp = eigen(Symmetric(ws.F_p[active_p, active_p]), Symmetric(ws.S[active_p, active_p]))
        
        # --- UPDATE ORBITALS ---
        # 1s
        c_1s_new = zeros(Float64, n); c_1s_new[active_s] = evecs_fs[:, 1]
        c_1s_new ./= sqrt(dot(c_1s_new, ws.S * c_1s_new))
        orbitals[1].coeffs = MIXING * orbitals[1].coeffs + (1 - MIXING) * c_1s_new
        orbitals[1].coeffs ./= sqrt(dot(orbitals[1].coeffs, ws.S * orbitals[1].coeffs))
        orbitals[1].energy = evals_fs[1]

        # 2s
        c_2s_new = zeros(Float64, n); c_2s_new[active_s] = evecs_fs[:, 2]
        c_2s_new ./= sqrt(dot(c_2s_new, ws.S * c_2s_new))
        orbitals[2].coeffs = MIXING * orbitals[2].coeffs + (1 - MIXING) * c_2s_new
        orbitals[2].coeffs ./= sqrt(dot(orbitals[2].coeffs, ws.S * orbitals[2].coeffs))
        orbitals[2].energy = evals_fs[2]

        # 2p
        c_2p_new = zeros(Float64, n); c_2p_new[active_p] = evecs_fp[:, 1]
        c_2p_new ./= sqrt(dot(c_2p_new, ws.S * c_2p_new))
        orbitals[3].coeffs = MIXING * orbitals[3].coeffs + (1 - MIXING) * c_2p_new
        orbitals[3].coeffs ./= sqrt(dot(orbitals[3].coeffs, ws.S * orbitals[3].coeffs))
        orbitals[3].energy = evals_fp[1]
        
        # 3s
        c_3s_new = zeros(Float64, n); c_3s_new[active_s] = evecs_fs[:, 3]
        c_3s_new ./= sqrt(dot(c_3s_new, ws.S * c_3s_new))
        orbitals[4].coeffs = MIXING * orbitals[4].coeffs + (1 - MIXING) * c_3s_new
        orbitals[4].coeffs ./= sqrt(dot(orbitals[4].coeffs, ws.S * orbitals[4].coeffs))
        orbitals[4].energy = evals_fs[3]

        # 3p
        c_3p_new = zeros(Float64, n); c_3p_new[active_p] = evecs_fp[:, 2]
        c_3p_new ./= sqrt(dot(c_3p_new, ws.S * c_3p_new))
        orbitals[5].coeffs = MIXING * orbitals[5].coeffs + (1 - MIXING) * c_3p_new
        orbitals[5].coeffs ./= sqrt(dot(orbitals[5].coeffs, ws.S * orbitals[5].coeffs))
        orbitals[5].energy = evals_fp[2]
 
        # --- CALCULATE TOTAL ENERGY ---
        E_total = 0.0
        
        h_11 = dot(orbitals[1].coeffs, H_core_s * orbitals[1].coeffs)
        E_total += (orbitals[1].occ / 2.0) * (h_11 + orbitals[1].energy)
        
        h_22 = dot(orbitals[2].coeffs, H_core_s * orbitals[2].coeffs)
        E_total += (orbitals[2].occ / 2.0) * (h_22 + orbitals[2].energy)
        
        h_2p = dot(orbitals[3].coeffs, H_core_p * orbitals[3].coeffs)
        E_total += (orbitals[3].occ / 2.0) * (h_2p + orbitals[3].energy)
        
        h_33 = dot(orbitals[4].coeffs, H_core_s * orbitals[4].coeffs)
        E_total += (orbitals[4].occ / 2.0) * (h_33 + orbitals[4].energy)
        
        h_3p = dot(orbitals[5].coeffs, H_core_p * orbitals[5].coeffs)
        E_total += (orbitals[5].occ / 2.0) * (h_3p + orbitals[5].energy)
        
        delta = abs(E_total - E_old)
        t_iter = time() - t0
        
        @printf("%4d | %15.8f | %10.2e | %7.4fs\n", iter, E_total, delta, t_iter)
        
        if delta < 1e-9
            println("\nConverged!")
            println("Box size: $R_max a.u.")
            @printf("Final Energy: %.6f Ha\n", E_total)
            println("Reference   : ~ -526.817 Ha")
            println("===== END =====")
            break
        end
        
        E_old = E_total
    end
end

solve_argon(20.0)
