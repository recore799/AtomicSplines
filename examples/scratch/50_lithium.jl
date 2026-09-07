using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines

using LinearAlgebra
using Printf
using JLD2

function solve_lithium_rohf(R_max; verbose::Bool=true)
    println("=== Lithium ROHF Term-Dependent 2S (Z=3) ===")
    
    N_elems = 100
    Z = 3.0
    basis = generate_basis(R_max, N_elems, Val(7), γ=2.5)
    ws = init_scf_workspace(basis, Z)

    n = basis.num_splines
    active_s = 2:(n-1)  
    
    # 1s^2, 2s^1
    orbitals = [
        Orbital(1, 0, 2.0), # 1s (Closed Core)
        Orbital(2, 0, 1.0)  # 2s (Open Valence, w=1)
    ]
    
    H_core_s = ws.T + ws.V
    
    # --- Initial Guess ---
    evals_s, evecs_s = eigen(Symmetric(H_core_s[active_s, active_s]), ws.S[active_s, active_s])
    orbitals[1].coeffs = zeros(Float64, n); orbitals[1].coeffs[active_s] = evecs_s[:, 1]
    orbitals[2].coeffs = zeros(Float64, n); orbitals[2].coeffs[active_s] = evecs_s[:, 2]
    
    for orb in orbitals
        orb.coeffs ./= sqrt(dot(orb.coeffs, ws.S * orb.coeffs))
    end

    # Pre-allocate matrices
    J_1s = zeros(Float64, n, n)
    K_1s = zeros(Float64, n, n)
    
    J_2s = zeros(Float64, n, n)
    K_2s = zeros(Float64, n, n)
    
    F_1s = zeros(Float64, n, n)
    F_2s = zeros(Float64, n, n)

    E_old = 0.0
    MIXING = 0.3 
    
    println("Comenzando ciclo SCF ROHF para Litio...")
    if verbose
        @printf("%-4s | %-14s | %-10s | %-8s\n", "Iter", "E_total (Ha)", "Delta E", "Time (s)")
        println("-"^78)
    end

    for iter in 1:100
        t0 = time()

        # ==========================================
        # PHASE 1: Build Individual Shell Potentials
        # ==========================================
        
        # 1s Core spherical background (w = 2.0)
        build_specific_J_matrix!(ws, J_1s, orbitals[1], 2.0)
        assemble_K_matrix!(ws, ws.K_mats[0], 0, [orbitals[1]])
        K_1s .= ws.K_mats[0]
        
        # 2s Valence spherical background (w = 1.0)
        build_specific_J_matrix!(ws, J_2s, orbitals[2], 1.0)
        assemble_K_matrix!(ws, ws.K_mats[0], 0, [orbitals[2]])
        K_2s .= ws.K_mats[0]

        # ==========================================
        # PHASE 2: Assemble Shell-Dependent Fock Matrices
        # ==========================================
        
        # 1s Matrix: Feels the full atomic density. Self-interaction cancels exactly via Roothaan math.
        F_1s .= H_core_s .+ J_1s .+ J_2s .- K_1s .- K_2s
        
        # 2s Matrix: Feels the 1s core. It does NOT feel itself (w-1 = 0).
        F_2s .= H_core_s .+ J_1s .- K_1s
        
        # ==========================================
        # PHASE 3: Diagonalization & Gram-Schmidt
        # ==========================================
        
        evals_f1s, evecs_f1s = eigen(Symmetric(F_1s[active_s, active_s]), ws.S[active_s, active_s])
        evals_f2s, evecs_f2s = eigen(Symmetric(F_2s[active_s, active_s]), ws.S[active_s, active_s])
        
        # Extract new raw vectors
        c_1s_new = zeros(Float64, n); c_1s_new[active_s] = evecs_f1s[:, 1]
        
        # We must pull the second lowest energy state for 2s to maintain the radial node!
        c_2s_new = zeros(Float64, n); c_2s_new[active_s] = evecs_f2s[:, 2]
        
        # --- Gram-Schmidt Orthogonalization ---
        # The 1s core is rigid; we rigidly project it out of the sensitive 2s valence
        overlap_1s_2s = dot(c_1s_new, ws.S * c_2s_new)
        c_2s_new .-= overlap_1s_2s .* c_1s_new
        
        # Normalize
        c_1s_new ./= sqrt(dot(c_1s_new, ws.S * c_1s_new))
        c_2s_new ./= sqrt(dot(c_2s_new, ws.S * c_2s_new))
        
        # --- Update Orbitals ---
        orbitals[1].coeffs = MIXING * orbitals[1].coeffs + (1 - MIXING) * c_1s_new
        orbitals[2].coeffs = MIXING * orbitals[2].coeffs + (1 - MIXING) * c_2s_new
        
        for orb in orbitals
            orb.coeffs ./= sqrt(dot(orb.coeffs, ws.S * orb.coeffs))
        end

        orbitals[1].energy = evals_f1s[1]
        orbitals[2].energy = evals_f2s[2]
 
        # ==========================================
        # PHASE 4: Compute Total Energy
        # ==========================================
        E_total = 0.0

        for orb in orbitals
            if orb.occ > 0.0
                h_ii = dot(orb.coeffs, H_core_s * orb.coeffs)
                E_total += (orb.occ / 2.0) * (h_ii + orb.energy)
            end
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
            @printf("Energía final HF-t (^2S): %.6f Ha\n", E_total)

            filename = "lithium_rohf_results_R$(R_max).jld2"
            jldsave(filename;
                orbitals = orbitals,
                E_total = E_total,
                R_max = R_max,
                R_grid = ws.R,
                V_nuclear = ws.V,
                num_splines = n,
                active_s = active_s
            )
            println("Saved Term-Dependent ROHF data to $filename")
            println("===== END =====")
            break
        end
        
        E_old = E_total
    end
end

solve_lithium_rohf(30.0)
