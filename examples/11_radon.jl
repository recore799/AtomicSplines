using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using LinearAlgebra

using Printf
using JLD2

function build_active_density(orbitals, active_idx, target_l)
    dim = length(active_idx)
    D = zeros(Float64, dim, dim)
    for orb in orbitals
        if orb.l == target_l
            c_act = orb.coeffs[active_idx]
            # Outer product of coefficients scaled by occupation
            D .+= orb.occ .* (c_act * c_act') 
        end
    end
    return D
end

function solve_radon(R_max; verbose::Bool=true)
    println("=== Radon (Z=86) ===")
    
    time1 = time()
    N_elems = 400
    Z = 86.0
    basis = generate_basis(R_max, N_elems, Val(7), γ=2.5)
    ws = init_scf_workspace(basis, Z)

    n = basis.num_splines

    # Enforce origin boundary conditions: P(r) ~ r^(l+1)
    active_s = 2:(n-1)  # P(0) = 0
    active_p = 3:(n-1)  # P(0) = 0, P'(0) = 0
    active_d = 4:(n-1)  # P(0) = 0, P'(0) = 0, P''(0) = 0
    active_f = 5:(n-1)  # P(0) = 0, P'(0) = 0, P''(0) = 0, P'''(0) = 0
    
    if verbose
        println("\n--- Información de la Base ---")
        println("  Radio máx (R_max)  : $R_max a.u.")
        println("  Elementos          : $N_elems")
        println("  Splines totales (n): $n")
        println("  Funciones s activas: $(length(active_s))")
        println("  Funciones p activas: $(length(active_p))")
        println("  Funciones d activas: $(length(active_d))")
        println("  Funciones f activas: $(length(active_f))")
        println("------------------------------\n")
    end

    # Initialize Orbitals for Radon (15 distinct subshells)
    orbitals = [
        # s-orbitals (Indices 1 to 6)
        Orbital(1, 0, 2.0), Orbital(2, 0, 2.0), Orbital(3, 0, 2.0), 
        Orbital(4, 0, 2.0), Orbital(5, 0, 2.0), Orbital(6, 0, 2.0),
        
        # p-orbitals (Indices 7 to 11)
        Orbital(2, 1, 6.0), Orbital(3, 1, 6.0), Orbital(4, 1, 6.0), 
        Orbital(5, 1, 6.0), Orbital(6, 1, 6.0),
        
        # d-orbitals (Indices 12 to 14)
        Orbital(3, 2, 10.0), Orbital(4, 2, 10.0), Orbital(5, 2, 10.0),
        
        # f-orbitals (Index 15)
        Orbital(4, 3, 14.0)
    ]
    
    # --- Core Hamiltonians ---
    H_core_s = ws.T + ws.V
    H_core_p = ws.T + ws.V + ws.V2 
    H_core_d = ws.T + ws.V + 3.0 * ws.V2
    H_core_f = ws.T + ws.V + 6.0 * ws.V2  # l(l+1)/2 = 12/2 = 6
    
    # --- Initial Guess ---
    evals_s, evecs_s = eigen(Symmetric(H_core_s[active_s, active_s]), ws.S[active_s, active_s])
    evals_p, evecs_p = eigen(Symmetric(H_core_p[active_p, active_p]), ws.S[active_p, active_p])
    evals_d, evecs_d = eigen(Symmetric(H_core_d[active_d, active_d]), ws.S[active_d, active_d])
    evals_f, evecs_f = eigen(Symmetric(H_core_f[active_f, active_f]), ws.S[active_f, active_f])
    
    for (i, orb_idx) in enumerate(1:6) # 1s to 6s
        orbitals[orb_idx].coeffs = zeros(Float64, n)
        orbitals[orb_idx].coeffs[active_s] = evecs_s[:, i]
    end
    for (i, orb_idx) in enumerate(7:11) # 2p to 6p
        orbitals[orb_idx].coeffs = zeros(Float64, n)
        orbitals[orb_idx].coeffs[active_p] = evecs_p[:, i]
    end
    for (i, orb_idx) in enumerate(12:14) # 3d to 5d
        orbitals[orb_idx].coeffs = zeros(Float64, n)
        orbitals[orb_idx].coeffs[active_d] = evecs_d[:, i]
    end
    for (i, orb_idx) in enumerate(15:15) # 4f
        orbitals[orb_idx].coeffs = zeros(Float64, n)
        orbitals[orb_idx].coeffs[active_f] = evecs_f[:, i]
    end
    
    # Normalize initial guesses
    for orb in orbitals
        orb.coeffs ./= sqrt(dot(orb.coeffs, ws.S * orb.coeffs))
    end
    # --- SCF Loop ---
    E_old = 0.0
    MIXING = 0.3
    
    println("Comenzando ciclo SCF...")
    if verbose
        @printf("%-4s | %-14s | %-10s | %-8s\n", "Iter", "E_total (Ha)", "Delta E", "Time (s)")
        println("-"^78)
    end

    max_diis_history = 6  # 6 to 8 is usually the sweet spot for heavy atoms
    F_hist = [] 
    E_hist = []

    for iter in 1:1000
        t0 = time()

        build_total_J_matrix!(ws, orbitals)

        assemble_K_matrix!(ws, ws.K_mats[0], 0, orbitals)
        assemble_K_matrix!(ws, ws.K_mats[1], 1, orbitals)
        assemble_K_matrix!(ws, ws.K_mats[2], 2, orbitals)
        assemble_K_matrix!(ws, ws.K_mats[3], 3, orbitals) # Added for f
        
        # --- Build Fock Matrices ---
        ws.F_s .= H_core_s .+ ws.J .- ws.K_mats[0]
        ws.F_p .= H_core_p .+ ws.J .- ws.K_mats[1]
        ws.F_d .= H_core_d .+ ws.J .- ws.K_mats[2]
        ws.F_f .= H_core_f .+ ws.J .- ws.K_mats[3] # Added for f

        # --- Compute Density Matrices for Active Blocks ---
        D_s = build_active_density(orbitals, active_s, 0)
        D_p = build_active_density(orbitals, active_p, 1)
        D_d = build_active_density(orbitals, active_d, 2)
        D_f = build_active_density(orbitals, active_f, 3)

        # --- Compute Commutator Error Matrices (E = FDS - SDF) ---
        S_s, S_p = ws.S[active_s, active_s], ws.S[active_p, active_p]
        S_d, S_f = ws.S[active_d, active_d], ws.S[active_f, active_f]

        F_s_act, F_p_act = ws.F_s[active_s, active_s], ws.F_p[active_p, active_p]
        F_d_act, F_f_act = ws.F_d[active_d, active_d], ws.F_f[active_f, active_f]

        E_s = F_s_act * D_s * S_s - S_s * D_s * F_s_act
        E_p = F_p_act * D_p * S_p - S_p * D_p * F_p_act
        E_d = F_d_act * D_d * S_d - S_d * D_d * F_d_act
        E_f = F_f_act * D_f * S_f - S_f * D_f * F_f_act

        # --- Store History ---
        push!(F_hist, (copy(ws.F_s), copy(ws.F_p), copy(ws.F_d), copy(ws.F_f)))
        push!(E_hist, (E_s, E_p, E_d, E_f))

        if length(F_hist) > max_diis_history
            popfirst!(F_hist)
            popfirst!(E_hist)
        end

        num_hist = length(F_hist)
        
        # --- Build and Solve the Pulay Matrix ---
        if num_hist >= 2
            B = zeros(Float64, num_hist + 1, num_hist + 1)
            for i in 1:num_hist
                for j in i:num_hist
                    # The inner product is the sum of Frobenius products over all l-blocks
                    val = dot(E_hist[i][1], E_hist[j][1]) + 
                          dot(E_hist[i][2], E_hist[j][2]) + 
                          dot(E_hist[i][3], E_hist[j][3]) + 
                          dot(E_hist[i][4], E_hist[j][4])
                    B[i, j] = val
                    B[j, i] = val
                end
                B[i, end] = -1.0
                B[end, i] = -1.0
            end
            B[end, end] = 0.0

            rhs = zeros(Float64, num_hist + 1)
            rhs[end] = -1.0

            # Solve for the mixing coefficients c_i
            c_diis = B \ rhs

            # Extrapolate the new Fock matrices
            F_s_eff, F_p_eff = zeros(size(ws.F_s)), zeros(size(ws.F_p))
            F_d_eff, F_f_eff = zeros(size(ws.F_d)), zeros(size(ws.F_f))

            for i in 1:num_hist
                F_s_eff .+= c_diis[i] .* F_hist[i][1]
                F_p_eff .+= c_diis[i] .* F_hist[i][2]
                F_d_eff .+= c_diis[i] .* F_hist[i][3]
                F_f_eff .+= c_diis[i] .* F_hist[i][4]
            end
        else
            # Not enough history yet, use current iteration
            F_s_eff, F_p_eff = ws.F_s, ws.F_p
            F_d_eff, F_f_eff = ws.F_d, ws.F_f
        end

        # --- Diagonalize the Effective Fock Matrices ---
        evals_fs, evecs_fs = eigen(Symmetric(F_s_eff[active_s, active_s]), S_s)
        evals_fp, evecs_fp = eigen(Symmetric(F_p_eff[active_p, active_p]), S_p)
        evals_fd, evecs_fd = eigen(Symmetric(F_d_eff[active_d, active_d]), S_d)
        evals_ff, evecs_ff = eigen(Symmetric(F_f_eff[active_f, active_f]), S_f)

        # --- Update Orbitals (NO MIXING NEEDED) ---
        # With DIIS, you take the resulting eigenvectors exactly as they are.
        
        for (i, orb_idx) in enumerate(1:6) # s
            orbitals[orb_idx].coeffs .= 0.0
            orbitals[orb_idx].coeffs[active_s] = evecs_fs[:, i]
            orbitals[orb_idx].coeffs ./= sqrt(dot(orbitals[orb_idx].coeffs, ws.S * orbitals[orb_idx].coeffs))
            orbitals[orb_idx].energy = evals_fs[i]
        end

        for (i, orb_idx) in enumerate(7:11) # p
            orbitals[orb_idx].coeffs .= 0.0
            orbitals[orb_idx].coeffs[active_p] = evecs_fp[:, i]
            orbitals[orb_idx].coeffs ./= sqrt(dot(orbitals[orb_idx].coeffs, ws.S * orbitals[orb_idx].coeffs))
            orbitals[orb_idx].energy = evals_fp[i]
        end

        for (i, orb_idx) in enumerate(12:14) # d
            orbitals[orb_idx].coeffs .= 0.0
            orbitals[orb_idx].coeffs[active_d] = evecs_fd[:, i]
            orbitals[orb_idx].coeffs ./= sqrt(dot(orbitals[orb_idx].coeffs, ws.S * orbitals[orb_idx].coeffs))
            orbitals[orb_idx].energy = evals_fd[i]
        end

        for (i, orb_idx) in enumerate(15:15) # f
            orbitals[orb_idx].coeffs .= 0.0
            orbitals[orb_idx].coeffs[active_f] = evecs_ff[:, i]
            orbitals[orb_idx].coeffs ./= sqrt(dot(orbitals[orb_idx].coeffs, ws.S * orbitals[orb_idx].coeffs))
            orbitals[orb_idx].energy = evals_ff[i]
        end
        
        # --- Diagonalize Independent Blocks ---
        # evals_fs, evecs_fs = eigen(Symmetric(ws.F_s[active_s, active_s]), ws.S[active_s, active_s])
        # evals_fp, evecs_fp = eigen(Symmetric(ws.F_p[active_p, active_p]), ws.S[active_p, active_p])
        # evals_fd, evecs_fd = eigen(Symmetric(ws.F_d[active_d, active_d]), ws.S[active_d, active_d])
        # evals_ff, evecs_ff = eigen(Symmetric(ws.F_f[active_f, active_f]), ws.S[active_f, active_f])
        
        # # --- Update Orbitals ---
        
        # for (i, orb_idx) in enumerate(1:6)
        #     c_new = zeros(Float64, n); c_new[active_s] = evecs_fs[:, i]
        #     c_new ./= sqrt(dot(c_new, ws.S * c_new))
        #     orbitals[orb_idx].coeffs = MIXING * orbitals[orb_idx].coeffs + (1 - MIXING) * c_new
        #     orbitals[orb_idx].coeffs ./= sqrt(dot(orbitals[orb_idx].coeffs, ws.S * orbitals[orb_idx].coeffs))
        #     orbitals[orb_idx].energy = evals_fs[i]
        # end

        # for (i, orb_idx) in enumerate(7:11)
        #     c_new = zeros(Float64, n); c_new[active_p] = evecs_fp[:, i]
        #     c_new ./= sqrt(dot(c_new, ws.S * c_new))
        #     orbitals[orb_idx].coeffs = MIXING * orbitals[orb_idx].coeffs + (1 - MIXING) * c_new
        #     orbitals[orb_idx].coeffs ./= sqrt(dot(orbitals[orb_idx].coeffs, ws.S * orbitals[orb_idx].coeffs))
        #     orbitals[orb_idx].energy = evals_fp[i]
        # end

        # for (i, orb_idx) in enumerate(12:14)
        #     c_new = zeros(Float64, n); c_new[active_d] = evecs_fd[:, i]
        #     c_new ./= sqrt(dot(c_new, ws.S * c_new))
        #     orbitals[orb_idx].coeffs = MIXING * orbitals[orb_idx].coeffs + (1 - MIXING) * c_new
        #     orbitals[orb_idx].coeffs ./= sqrt(dot(orbitals[orb_idx].coeffs, ws.S * orbitals[orb_idx].coeffs))
        #     orbitals[orb_idx].energy = evals_fd[i]
        # end

        # for (i, orb_idx) in enumerate(15:15) # f-orbital update
        #     c_new = zeros(Float64, n); c_new[active_f] = evecs_ff[:, i]
        #     c_new ./= sqrt(dot(c_new, ws.S * c_new))
        #     orbitals[orb_idx].coeffs = MIXING * orbitals[orb_idx].coeffs + (1 - MIXING) * c_new
        #     orbitals[orb_idx].coeffs ./= sqrt(dot(orbitals[orb_idx].coeffs, ws.S * orbitals[orb_idx].coeffs))
        #     orbitals[orb_idx].energy = evals_ff[i]
        # end

        # --- Compute Total Energy ---
        E_total = 0.0
        
        for orb in orbitals
            if orb.l == 0
                h_ii = dot(orb.coeffs, H_core_s * orb.coeffs)
            elseif orb.l == 1
                h_ii = dot(orb.coeffs, H_core_p * orb.coeffs)
            elseif orb.l == 2
                h_ii = dot(orb.coeffs, H_core_d * orb.coeffs)
            elseif orb.l == 3
                h_ii = dot(orb.coeffs, H_core_f * orb.coeffs)
            end
            E_total += (orb.occ / 2.0) * (h_ii + orb.energy)
        end
        
        delta = abs(E_total - E_old)
        elapsed = time() - t0
        
        if verbose
            @printf("%-4d | %14.8f | %10.2e | %8.4f\n", iter, E_total, delta, elapsed)
        end
        
        if delta < 1e-9
            if verbose 
                println("-"^78)
                println("Converged in $iter iterations.")
            end
            time_total = time() - time1
            println("Radio de confinamiento: $R_max a.u.")
            @printf("Energía final: %.6f Ha\n", E_total)
            println("Límite Hartree-Fock: -21866.772 Ha") 
            @printf("Tiempo elapsado: %.4f s\n", time_total)
            println("===== END =====")

            # --- Save the results ---
            filename = "radon_results_R$(R_max).jld2"
            @save filename orbitals E_total R_max
            println("Saved orbitals and energy to $filename")
            
            # return orbitals, E_total # It's good practice to return them too
            break
        end
        
        E_old = E_total
    end
end


solve_radon(20.0)
