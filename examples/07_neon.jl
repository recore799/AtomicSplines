using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using LinearAlgebra

using Printf


function solve_neon(R_max; verbose::Bool=true)
    println("=== Neon (Z=10) ===")
    
    N_elems = 200
    Z = 10.0
    basis = generate_basis(R_max, N_elems, Val(7), γ=2.5)
    ws = init_scf_workspace(basis, Z)

    n = basis.num_splines

    # Enforce origin boundary conditions: P(r) ~ r^(l+1)
    active_s = 2:(n-1)  # s-orbitals (l=0): Drop B1 to enforce P(0) = 0
    active_p = 3:(n-1)  # p-orbitals (l=1): Drop B1 and B2 to enforce P(0)=0, P'(0)=0
    
    if verbose
        println("\n--- Información de la Base ---")
        println("  Radio máx (R_max)  : $R_max a.u.")
        println("  Elementos          : $N_elems")
        println("  Splines totales (n): $n")
        println("  Funciones s activas: $(length(active_s))")
        println("  Funciones p activas: $(length(active_p))")
        println("------------------------------\n")
    end

    # Initialize Orbitals (n, l, occupancy)
    orbitals = [
        Orbital(1, 0, 2.0), # 1s
        Orbital(2, 0, 2.0), # 2s
        Orbital(2, 1, 6.0)  # 2p
    ]
    
    # --- Core Hamiltonians ---
    H_core_s = ws.T + ws.V
    H_core_p = ws.T + ws.V + ws.V2 
    
    # --- Initial Guess ---
    evals_s, evecs_s = eigen(Symmetric(H_core_s[active_s, active_s]), ws.S[active_s, active_s])
    orbitals[1].coeffs = zeros(Float64, n); orbitals[1].coeffs[active_s] = evecs_s[:, 1]
    orbitals[2].coeffs = zeros(Float64, n); orbitals[2].coeffs[active_s] = evecs_s[:, 2]
    
    evals_p, evecs_p = eigen(Symmetric(H_core_p[active_p, active_p]), ws.S[active_p, active_p])
    orbitals[3].coeffs = zeros(Float64, n); orbitals[3].coeffs[active_p] = evecs_p[:, 1]
    
    # Normalize initial guesses
    for orb in orbitals
        orb.coeffs ./= sqrt(dot(orb.coeffs, ws.S * orb.coeffs))
    end

    # --- SCF Loop ---
    E_old = 0.0
    MIXING = 0.0 # Dampening factor to aid convergence
    
    println("Comenzando ciclo SCF...")
    if verbose
        # Added Time (s) to the header and expanded the separator line
        @printf("%-4s | %-14s | %-10s | %-8s\n", 
                "Iter", "E_total (Ha)", "Delta E", "Time (s)")
        println("-"^78)
    end

    for iter in 1:60
        t0 = time()

        build_total_J_matrix!(ws, orbitals)

        assemble_K_matrix!(ws, ws.K_mats[0] , 0, orbitals)
        assemble_K_matrix!(ws, ws.K_mats[1] , 1, orbitals)
        
        # --- Build Fock Matrices ---
        ws.F_s .= H_core_s .+ ws.J .- ws.K_mats[0]
        ws.F_p .= H_core_p .+ ws.J .- ws.K_mats[1]
        
        # --- Diagonalize Independent Blocks ---
        evals_fs, evecs_fs = eigen(Symmetric(ws.F_s[active_s, active_s]), ws.S[active_s, active_s])
        evals_fp, evecs_fp = eigen(Symmetric(ws.F_p[active_p, active_p]), ws.S[active_p, active_p])
        
        # --- Update Orbitals (With Mixing) ---
        # 1s
        c_1s_new = zeros(Float64, n); c_1s_new[active_s] = evecs_fs[:, 1]
        c_1s_new ./= sqrt(dot(c_1s_new, ws.S * c_1s_new))
        # if dot(orbitals[1].coeffs, ws.S * c_1s_new) < 0
        #     c_1s_new .= -c_1s_new 
        # end
        orbitals[1].coeffs = MIXING * orbitals[1].coeffs + (1 - MIXING) * c_1s_new
        orbitals[1].coeffs ./= sqrt(dot(orbitals[1].coeffs, ws.S * orbitals[1].coeffs))
        orbitals[1].energy = evals_fs[1]

        # 2s
        c_2s_new = zeros(Float64, n); c_2s_new[active_s] = evecs_fs[:, 2]
        c_2s_new ./= sqrt(dot(c_2s_new, ws.S * c_2s_new))
        # if dot(orbitals[2].coeffs, ws.S * c_2s_new) < 0
        #     c_2s_new .= -c_2s_new
        # end
        orbitals[2].coeffs = MIXING * orbitals[2].coeffs + (1 - MIXING) * c_2s_new
        orbitals[2].coeffs ./= sqrt(dot(orbitals[2].coeffs, ws.S * orbitals[2].coeffs))
        orbitals[2].energy = evals_fs[2]

        # 2p
        c_2p_new = zeros(Float64, n); c_2p_new[active_p] = evecs_fp[:, 1]
        c_2p_new ./= sqrt(dot(c_2p_new, ws.S * c_2p_new))
        # if dot(orbitals[3].coeffs, ws.S * c_2p_new) < 0
        #     c_2p_new .= -c_2p_new
        # end
        orbitals[3].coeffs = MIXING * orbitals[3].coeffs + (1 - MIXING) * c_2p_new
        orbitals[3].coeffs ./= sqrt(dot(orbitals[3].coeffs, ws.S * orbitals[3].coeffs))
        orbitals[3].energy = evals_fp[1]
 
        # --- Compute Total Energy ---
        E_total = 0.0


        h_11 = dot(orbitals[1].coeffs, H_core_s * orbitals[1].coeffs)
        E_total += (orbitals[1].occ / 2.0) * (h_11 + orbitals[1].energy)
        
        h_22 = dot(orbitals[2].coeffs, H_core_s * orbitals[2].coeffs)
        E_total += (orbitals[2].occ / 2.0) * (h_22 + orbitals[2].energy)
        
        h_pp = dot(orbitals[3].coeffs, H_core_p * orbitals[3].coeffs)
        E_total += (orbitals[3].occ / 2.0) * (h_pp + orbitals[3].energy)
        
        delta = abs(E_total - E_old)
        
        # Calculate iteration time
        elapsed = time() - t0
        
        if verbose
            @printf("%-4d | %14.8f | %10.2e | %8.4f\n", 
                    iter, E_total, delta, elapsed)
        end
        
        if delta < 1e-9
            if verbose 
                println("-"^78)
                println("Converged in $iter iterations.")
            end
            println("Radio de confinamiento: $R_max a.u.")
            @printf("Energía final: %.6f Ha\n", E_total)
            println("===== END =====")
            break
        end
        
        E_old = E_total
    end
end

solve_neon(1.0)
