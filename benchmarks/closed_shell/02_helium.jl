using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using LinearAlgebra
using Printf


function solve_helium(R_max; verbose::Bool=false)
    println("=== Helium (Z=2) ===")
    
    N_elems = 30
    Z = 2.0 # Changed from 10.0 to 2.0 for Helium
    basis = generate_basis(R_max, N_elems, Val(7), γ=2.5)
    ws = init_scf_workspace(basis, Z)

    n = basis.num_splines

    # Enforce origin boundary conditions: P(r) ~ r^(l+1)
    active_s = 2:(n-1)  # s-orbitals (l=0): Drop B1 to enforce P(0) = 0

   
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
        
        # --- Build Fock Matrix ---
        F_s = ws.F_mats[0]
        F_s .= H_core_s .+ ws.J .- ws.K_mats[0]
        
        # --- Diagonalize Independent Blocks ---
        evals_fs, evecs_fs = eigen(Symmetric(F_s[active_s, active_s]), ws.S[active_s, active_s])
        
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
            println("Radio de confinamiento: $R_max a.u.")
            @printf("Energía final: %.6f Ha\n", E_total)
            println("===== END =====")
            
            break
        end
        
        E_old = E_total
    end
end

solve_helium(10.0, verbose=true)

