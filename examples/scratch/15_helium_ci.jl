using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using LinearAlgebra
using Printf
using JLD2


function solve_helium(R_max; verbose::Bool=false)
    println("=== Helium (Z=2) ===")
    
    N_elems = 500
    Z = 2.0 # Changed from 10.0 to 2.0 for Helium
    basis = generate_basis(R_max, N_elems, Val(7), γ=2.5)
    ws = init_scf_workspace(basis, Z)

    n = basis.num_splines

    # Enforce origin boundary conditions: P(r) ~ r^(l+1)
    active_s = 2:(n-1)
    active_p = 3:(n-1)
    active_d = 4:(n-1)
    active_f = 5:(n-1)
    active_g = 6:(n-1)
    active_i = 7:(n-1)

   
    if verbose
        println("\n--- Información de la Base ---")
        println("  Radio máx (R_max)  : $R_max a.u.")
        println("  Elementos          : $N_elems")
        println("  Splines totales (n): $n")
        println("  Funciones s activas: $(length(active_s))")
        println("------------------------------\n")
    end

    # Initialize Orbitals (n, l, occupancy)
    # Helium only has a fully occupied 1s orbital
    orbitals = [
        Orbital(1, 0, 2.0)  # 1s
    ]
    
    # --- Core Hamiltonians ---
    H_core_s = ws.T + ws.V

    H_core_p = ws.T + ws.V + ws.V2 
    H_core_d = ws.T + ws.V + 3.0 * ws.V2
    H_core_f = ws.T + ws.V + 6.0 * ws.V2
    H_core_g = ws.T + ws.V + 10.0 * ws.V2
    H_core_i = ws.T + ws.V + 15.0 * ws.V2

    # --- Initial Guess ---
    evals_s, evecs_s = eigen(Symmetric(H_core_s[active_s, active_s]), ws.S[active_s, active_s])
    orbitals[1].coeffs = zeros(Float64, n); orbitals[1].coeffs[active_s] = evecs_s[:, 1]
    
    orbitals[1].coeffs ./= sqrt(dot(orbitals[1].coeffs, ws.S * orbitals[1].coeffs))

    # --- SCF Loop ---
    E_old = 0.0
    MIXING = 0.0 # Dampening factor to aid convergence
    
    println("Comenzando ciclo SCF...")
    if verbose
        @printf("%-4s | %-14s | %-10s | %-8s\n", 
                "Iter", "E_total (Ha)", "Delta E", "Time (s)")
        println("-"^78)
    end

    for iter in 1:60
        t0 = time()

        build_total_J_matrix!(ws, orbitals)

        # Only need l=0 exchange matrix for Helium
        assemble_K_matrix!(ws, ws.K_mats[0] , 0, orbitals)
        
        # --- Build Fock Matrices ---
        ws.F_s .= H_core_s .+ ws.J .- ws.K_mats[0]
        
        # --- Diagonalize Independent Blocks ---
        evals_fs, evecs_fs = eigen(Symmetric(ws.F_s[active_s, active_s]), ws.S[active_s, active_s])
        
        # --- Update Orbitals (With Mixing) ---
        # 1s
        c_1s_new = zeros(Float64, n); c_1s_new[active_s] = evecs_fs[:, 1]
        c_1s_new ./= sqrt(dot(c_1s_new, ws.S * c_1s_new))
        
        orbitals[1].coeffs = MIXING * orbitals[1].coeffs + (1 - MIXING) * c_1s_new
        orbitals[1].coeffs ./= sqrt(dot(orbitals[1].coeffs, ws.S * orbitals[1].coeffs))
        orbitals[1].energy = evals_fs[1]
 
        # --- Compute Total Energy ---
        E_total = 0.0

        h_11 = dot(orbitals[1].coeffs, H_core_s * orbitals[1].coeffs)
        E_total += (orbitals[1].occ / 2.0) * (h_11 + orbitals[1].energy)
        
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
            # ==========================================
            # POST-SCF: GENERATE VIRTUAL BASIS FOR CI
            # ==========================================
            println("Generando orbitales virtuales p, d y f...")
            
            # 1. Build Exchange matrices
            assemble_K_matrix!(ws, ws.K_mats[1], 1, orbitals)
            assemble_K_matrix!(ws, ws.K_mats[2], 2, orbitals)
            assemble_K_matrix!(ws, ws.K_mats[3], 3, orbitals) 
            assemble_K_matrix!(ws, ws.K_mats[4], 4, orbitals) 
            assemble_K_matrix!(ws, ws.K_mats[5], 5, orbitals) 
            
            # 2. Build the exact Fock matrices
            ws.F_p .= H_core_p .+ ws.J .- ws.K_mats[1]
            ws.F_d .= H_core_d .+ ws.J .- ws.K_mats[2]
            ws.F_f .= H_core_f .+ ws.J .- ws.K_mats[3]
            ws.F_f .= H_core_f .+ ws.J .- ws.K_mats[4]
            ws.F_f .= H_core_f .+ ws.J .- ws.K_mats[5]
            
            # 3. Diagonalize to get the virtual eigenvectors
            evals_fp, evecs_fp = eigen(Symmetric(ws.F_p[active_p, active_p]), ws.S[active_p, active_p])
            evals_fd, evecs_fd = eigen(Symmetric(ws.F_d[active_d, active_d]), ws.S[active_d, active_d])
            evals_ff, evecs_ff = eigen(Symmetric(ws.F_f[active_f, active_f]), ws.S[active_f, active_f])
            evals_fg, evecs_fg = eigen(Symmetric(ws.F_g[active_g, active_g]), ws.S[active_g, active_g])
            evals_fi, evecs_fi = eigen(Symmetric(ws.F_i[active_i, active_i]), ws.S[active_i, active_i])
            
            # 4. Extract them
            N_VIRTUALS = 60
            # Keep s and p starting from the bottom
            virt_s = extract_virtuals(evals_fs, evecs_fs, 0, 1, N_VIRTUALS, active_s, n, ws)
            virt_p = extract_virtuals(evals_fp, evecs_fp, 1, 0, N_VIRTUALS, active_p, n, ws)

            virt_d = extract_virtuals(evals_fd, evecs_fd, 2, 0, N_VIRTUALS, active_d, n, ws, offset=2)

            virt_f = extract_virtuals(evals_ff, evecs_ff, 3, 0, N_VIRTUALS, active_f, n, ws, offset=15)
            virt_g = extract_virtuals(evals_fg, evecs_fg, 4, 0, N_VIRTUALS, active_g, n, ws, offset=25)
            virt_i = extract_virtuals(evals_fi, evecs_fi, 5, 0, N_VIRTUALS, active_i, n, ws, offset=35)


            # Combine them all
            all_virtuals = vcat(virt_s, virt_p, virt_d, virt_f, virt_g, virt_i)

            println("Total de orbitales virtuales extraídos: $(length(all_virtuals))")
            
            # Save orbitals and all_virtuals to a .jld2 file here
            filename = "helium_results_R$(R_max).jld2"
            @save filename orbitals all_virtuals E_total R_max
            println("Saved orbitals and energy to $filename")
            # ==========================================
            
            println("Radio de confinamiento: $R_max a.u.")
            @printf("Energía final: %.6f Ha\n", E_total)
            println("===== END =====")
            break
        end
        
        E_old = E_total
    end
end

solve_helium(20.0, verbose=true)

