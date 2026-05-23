using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using LinearAlgebra
using Printf
using JLD2


function solve_germanium_rohf(R_max; verbose::Bool=true)
    println("=== Germanium ROHF Term-Dependent 3P (Z=32) ===")
    
    N_elems = 300
    Z = 32.0

    ws = cached_init_scf_workspace(R_max, N_elems, Val(8), Z; γ=3.0, calc_R_matrices=true)
    basis = ws.basis

    n = basis.num_splines
    active_s = 2:(n-1)  
    active_p = 3:(n-1)  
    active_d = 4:(n-1)
    
    orbitals = [
        Orbital(1, 0, 2.0),  # [1] 1s (Closed)
        Orbital(2, 0, 2.0),  # [2] 2s (Closed)
        Orbital(2, 1, 6.0),  # [3] 2p (Closed)
        Orbital(3, 0, 2.0),  # [4] 3s (Closed)
        Orbital(3, 1, 6.0),  # [5] 3p (Closed)
        Orbital(3, 2, 10.0), # [6] 3d (Closed Core)
        Orbital(4, 0, 2.0),  # [7] 4s (Closed)
        Orbital(4, 1, 2.0)   # [8] 4p (Open Valence)
    ]
    
    H_core_s = ws.T + ws.V
    H_core_p = ws.T + ws.V + ws.R_inv2 
    H_core_d = ws.T + ws.V + 3.0 * ws.R_inv2 
    
    # --- Initial Guess ---
    evals_s, evecs_s = eigen(Symmetric(H_core_s[active_s, active_s]), ws.S[active_s, active_s])
    orbitals[1].coeffs = zeros(Float64, n); orbitals[1].coeffs[active_s] = evecs_s[:, 1]
    orbitals[2].coeffs = zeros(Float64, n); orbitals[2].coeffs[active_s] = evecs_s[:, 2]
    orbitals[4].coeffs = zeros(Float64, n); orbitals[4].coeffs[active_s] = evecs_s[:, 3]
    orbitals[7].coeffs = zeros(Float64, n); orbitals[7].coeffs[active_s] = evecs_s[:, 4]
    
    evals_p, evecs_p = eigen(Symmetric(H_core_p[active_p, active_p]), ws.S[active_p, active_p])
    orbitals[3].coeffs = zeros(Float64, n); orbitals[3].coeffs[active_p] = evecs_p[:, 1]
    orbitals[5].coeffs = zeros(Float64, n); orbitals[5].coeffs[active_p] = evecs_p[:, 2]
    orbitals[8].coeffs = zeros(Float64, n); orbitals[8].coeffs[active_p] = evecs_p[:, 3]

    evals_d, evecs_d = eigen(Symmetric(H_core_d[active_d, active_d]), ws.S[active_d, active_d])
    orbitals[6].coeffs = zeros(Float64, n); orbitals[6].coeffs[active_d] = evecs_d[:, 1]
    
    for orb in orbitals
        orb.coeffs ./= sqrt(dot(orb.coeffs, ws.S * orb.coeffs))
    end

    # Pre-allocate temporary matrices
    J_core = zeros(Float64, n, n); K_core_s = zeros(Float64, n, n)
    K_core_p = zeros(Float64, n, n); K_core_d = zeros(Float64, n, n)
    J_4p_spherical = zeros(Float64, n, n); K_4p_on_s = zeros(Float64, n, n)
    K_4p_on_p = zeros(Float64, n, n); K_4p_on_d = zeros(Float64, n, n)
    J_4p_intra = zeros(Float64, n, n); K_4p_intra = zeros(Float64, n, n)
    
    F_core_p = zeros(Float64, n, n); F_3d = zeros(Float64, n, n); F_4p = zeros(Float64, n, n)

    E_old = 0.0
    MIXING = 0.0
    
    println("Comenzando ciclo SCF ROHF...")
    if verbose
        @printf("%-4s | %-14s | %-10s | %-8s\n", "Iter", "E_total (Ha)", "Delta E", "Time (s)")
        println("-"^78)
    end

    for iter in 1:1000
        t0 = time()

        # --- Phase 1: Closed Core Matrix Assembly ---
        core_orbs = orbitals[1:7]
        build_total_J_matrix!(ws, core_orbs); J_core .= ws.J
        assemble_K_matrix!(ws, ws.K_mats[0], 0, core_orbs); K_core_s .= ws.K_mats[0]
        assemble_K_matrix!(ws, ws.K_mats[1], 1, core_orbs); K_core_p .= ws.K_mats[1]
        assemble_K_matrix!(ws, ws.K_mats[2], 2, core_orbs); K_core_d .= ws.K_mats[2]

        # --- Phase 2: Open 4p Shell Contributions ---
        build_specific_J_matrix!(ws, J_4p_spherical, orbitals[8], 2.0)
        assemble_K_matrix!(ws, ws.K_mats[0], 0, [orbitals[8]]); K_4p_on_s .= ws.K_mats[0]
        assemble_K_matrix!(ws, ws.K_mats[1], 1, [orbitals[8]]); K_4p_on_p .= ws.K_mats[1]
        assemble_K_matrix!(ws, ws.K_mats[2], 2, [orbitals[8]]); K_4p_on_d .= ws.K_mats[2]

        build_specific_J_matrix!(ws, J_4p_intra, orbitals[8], 1.0)
        build_specific_K_matrix!(ws, K_4p_intra, orbitals[8], 2, 5.0 / 25.0)

        # --- Phase 3: Raw Fock Matrices ---
        F_s = ws.F_mats[0]
        F_s .= H_core_s .+ J_core .+ J_4p_spherical .- K_core_s .- K_4p_on_s
        F_core_p .= H_core_p .+ J_core .+ J_4p_spherical .- K_core_p .- K_4p_on_p
        F_3d .= H_core_d .+ J_core .+ J_4p_spherical .- K_core_d .- K_4p_on_d
        F_4p .= H_core_p .+ J_core .+ J_4p_intra .- K_core_p .- K_4p_intra

        # ==========================================
        # PHASE 4: Scheduled Dynamic Level Shifting
        # ==========================================
        # Linearly decay the 3.0 Ha barrier over the first 20 iterations.
        level_shift = 3.0

        if level_shift > 0.0
            # Construct purely spatial geometric projectors
            P_s_spatial = zeros(Float64, n, n)
            for i in [1, 2, 4, 7] # 1s, 2s, 3s, 4s
                P_s_spatial .+= orbitals[i].coeffs * orbitals[i].coeffs'
            end
            
            P_p_spatial = zeros(Float64, n, n)
            for i in [3, 5, 8] # 2p, 3p, 4p
                P_p_spatial .+= orbitals[i].coeffs * orbitals[i].coeffs'
            end
            
            P_d_spatial = orbitals[6].coeffs * orbitals[6].coeffs' # 3d

            # Compute covariant projection matrices
            S_Ps_S = ws.S * P_s_spatial * ws.S
            S_Pp_S = ws.S * P_p_spatial * ws.S
            S_Pd_S = ws.S * P_d_spatial * ws.S

            # Apply shift to the virtual space
            F_s      .+= level_shift .* (ws.S .- S_Ps_S)
            F_core_p .+= level_shift .* (ws.S .- S_Pp_S)
            F_3d     .+= level_shift .* (ws.S .- S_Pd_S)
            F_4p     .+= level_shift .* (ws.S .- S_Pp_S)
            
            if verbose && iter <= 21
                @printf("  [Level Shift Active: b = %.3f Ha]\n", level_shift)
            end
        end

        # --- Phase 5: Diagonalization & Gram-Schmidt ---
        evals_fs, evecs_fs = eigen(Symmetric(F_s[active_s, active_s]), ws.S[active_s, active_s])
        evals_fcp, evecs_fcp = eigen(Symmetric(F_core_p[active_p, active_p]), ws.S[active_p, active_p])
        evals_f3d, evecs_f3d = eigen(Symmetric(F_3d[active_d, active_d]), ws.S[active_d, active_d])
        evals_f4p, evecs_f4p = eigen(Symmetric(F_4p[active_p, active_p]), ws.S[active_p, active_p])
        
        c_1s_new = zeros(Float64, n); c_1s_new[active_s] = evecs_fs[:, 1]
        c_2s_new = zeros(Float64, n); c_2s_new[active_s] = evecs_fs[:, 2]
        c_3s_new = zeros(Float64, n); c_3s_new[active_s] = evecs_fs[:, 3]
        c_4s_new = zeros(Float64, n); c_4s_new[active_s] = evecs_fs[:, 4]
        
        c_2p_new = zeros(Float64, n); c_2p_new[active_p] = evecs_fcp[:, 1]
        c_3p_new = zeros(Float64, n); c_3p_new[active_p] = evecs_fcp[:, 2]
        c_4p_new = zeros(Float64, n); c_4p_new[active_p] = evecs_f4p[:, 3]
        c_3d_new = zeros(Float64, n); c_3d_new[active_d] = evecs_f3d[:, 1]
        
        # Cascaded Gram-Schmidt Orthogonalization
        c_2s_new .-= dot(c_1s_new, ws.S * c_2s_new) .* c_1s_new
        c_3s_new .-= dot(c_1s_new, ws.S * c_3s_new) .* c_1s_new
        c_3s_new .-= dot(c_2s_new, ws.S * c_3s_new) .* c_2s_new
        c_4s_new .-= dot(c_1s_new, ws.S * c_4s_new) .* c_1s_new
        c_4s_new .-= dot(c_2s_new, ws.S * c_4s_new) .* c_2s_new
        c_4s_new .-= dot(c_3s_new, ws.S * c_4s_new) .* c_3s_new
        
        c_3p_new .-= dot(c_2p_new, ws.S * c_3p_new) .* c_2p_new
        c_4p_new .-= dot(c_2p_new, ws.S * c_4p_new) .* c_2p_new
        c_4p_new .-= dot(c_3p_new, ws.S * c_4p_new) .* c_3p_new
        
        c_1s_new ./= sqrt(dot(c_1s_new, ws.S * c_1s_new)); c_2s_new ./= sqrt(dot(c_2s_new, ws.S * c_2s_new))
        c_3s_new ./= sqrt(dot(c_3s_new, ws.S * c_3s_new)); c_4s_new ./= sqrt(dot(c_4s_new, ws.S * c_4s_new))
        c_2p_new ./= sqrt(dot(c_2p_new, ws.S * c_2p_new)); c_3p_new ./= sqrt(dot(c_3p_new, ws.S * c_3p_new))
        c_4p_new ./= sqrt(dot(c_4p_new, ws.S * c_4p_new)); c_3d_new ./= sqrt(dot(c_3d_new, ws.S * c_3d_new))
        
        # Linear Damping
        orbitals[1].coeffs = MIXING * orbitals[1].coeffs + (1 - MIXING) * c_1s_new
        orbitals[2].coeffs = MIXING * orbitals[2].coeffs + (1 - MIXING) * c_2s_new
        orbitals[3].coeffs = MIXING * orbitals[3].coeffs + (1 - MIXING) * c_2p_new
        orbitals[4].coeffs = MIXING * orbitals[4].coeffs + (1 - MIXING) * c_3s_new
        orbitals[5].coeffs = MIXING * orbitals[5].coeffs + (1 - MIXING) * c_3p_new
        orbitals[6].coeffs = MIXING * orbitals[6].coeffs + (1 - MIXING) * c_3d_new
        orbitals[7].coeffs = MIXING * orbitals[7].coeffs + (1 - MIXING) * c_4s_new
        orbitals[8].coeffs = MIXING * orbitals[8].coeffs + (1 - MIXING) * c_4p_new
        
        for orb in orbitals; orb.coeffs ./= sqrt(dot(orb.coeffs, ws.S * orb.coeffs)); end

        orbitals[1].energy = evals_fs[1]; orbitals[2].energy = evals_fs[2]; orbitals[3].energy = evals_fcp[1]
        orbitals[4].energy = evals_fs[3]; orbitals[5].energy = evals_fcp[2]; orbitals[6].energy = evals_f3d[1]
        orbitals[7].energy = evals_fs[4]; orbitals[8].energy = evals_f4p[3] 
 
        # --- Phase 6: Compute Total Energy ---
        E_total = 0.0
        for orb in orbitals
            if orb.occ > 0.0
                h_core = orb.l == 0 ? H_core_s : (orb.l == 1 ? H_core_p : H_core_d)
                E_total += (orb.occ / 2.0) * (dot(orb.coeffs, h_core * orb.coeffs) + orb.energy)
            end
        end
        
        delta = abs(E_total - E_old)
        elapsed = time() - t0
        
        if verbose
            @printf("%-4d | %14.8f | %10.2e | %8.4f\n", iter, E_total, delta, elapsed)
        end
        
        # Pure Energy Stabilization Metric
        if delta < 1e-9 
            if verbose 
                println("-"^78)
                println("Converged purely on energy metric in $iter iterations.")
            end
            @printf("Energía final HF-t (^3P): %.6f Ha\n", E_total)

            dense_grid = exp.(range(log(1e-8), log(R_max), length=10000))
            V_eff = compute_effective_central_potential(ws, orbitals, dense_grid, Z)
            P_4p = evaluate_orbital(ws.basis, orbitals[8].coeffs, dense_grid)

            filename = "germanium_rohf_results_R$(R_max).jld2"
            
            # Corrected export assignment mapping the geometric coordinates
            jldsave(filename;
                orbitals = orbitals, E_total = E_total, R_max = R_max,
                R_grid = dense_grid, # Exclusively continuous mapping
                V_nuclear = ws.V, num_splines = n, active_s = active_s,
                active_p = active_p, active_d = active_d,
                V_eff = V_eff, P_4p = P_4p
            )
            println("Saved Term-Dependent ROHF data and SO coupling prerequisites to $filename")
            println("===== END =====")
            break
        end
        
        E_old = E_total
    end
end

solve_germanium_rohf(30.0)
