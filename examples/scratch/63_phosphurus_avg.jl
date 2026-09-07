using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using LinearAlgebra
using Printf
using JLD2

function solve_phosphorus_rohf_CA(R_max; verbose::Bool=true)
    println("=== Phosphorus ROHF Configuration Average (Z=15) ===")
    
    N_elems = 100
    Z = 15.0
   
    ws = cached_init_scf_workspace(R_max, N_elems, Val(7), Z; γ=2.5, calc_R_matrices=true)
    
    basis = ws.basis
    n = basis.num_splines
    active_s = 2:(n-1)  
    active_p = 3:(n-1)  
   
    # 1s^2, 2s^2, 2p^6, 3s^2, 3p^3
    orbitals = [
        Orbital(1, 0, 2.0), # 1s (Closed)
        Orbital(2, 0, 2.0), # 2s (Closed)
        Orbital(2, 1, 6.0), # 2p (Closed Core)
        Orbital(3, 0, 2.0), # 3s (Closed)
        Orbital(3, 1, 3.0)  # 3p (Open Valence, fractional w=3)
    ]
    
    H_core_s = ws.T + ws.V
    H_core_p = ws.T + ws.V + ws.R_inv2 
    
    # --- Initial Guess ---
    evals_s, evecs_s = eigen(Symmetric(H_core_s[active_s, active_s]), ws.S[active_s, active_s])
    orbitals[1].coeffs = zeros(Float64, n);
    orbitals[1].coeffs[active_s] = evecs_s[:, 1]
    orbitals[2].coeffs = zeros(Float64, n);
    orbitals[2].coeffs[active_s] = evecs_s[:, 2]
    orbitals[4].coeffs = zeros(Float64, n);
    orbitals[4].coeffs[active_s] = evecs_s[:, 3]
    
    evals_p, evecs_p = eigen(Symmetric(H_core_p[active_p, active_p]), ws.S[active_p, active_p])
    orbitals[3].coeffs = zeros(Float64, n);
    orbitals[3].coeffs[active_p] = evecs_p[:, 1]
    orbitals[5].coeffs = zeros(Float64, n);
    orbitals[5].coeffs[active_p] = evecs_p[:, 2]
    
    for orb in orbitals
        orb.coeffs ./= sqrt(dot(orb.coeffs, ws.S * orb.coeffs))
    end

    # Pre-allocate specific temporary matrices
    J_core         = zeros(Float64, n, n)
    K_core_s       = zeros(Float64, n, n)
    K_core_p       = zeros(Float64, n, n)
    
    J_3p_spherical = zeros(Float64, n, n)
    K_3p_on_s      = zeros(Float64, n, n)
    K_3p_on_p      = zeros(Float64, n, n)
    
    J_3p_intra     = zeros(Float64, n, n)
    K_3p_intra     = zeros(Float64, n, n)
    
    F_2p = zeros(Float64, n, n)
    F_3p = zeros(Float64, n, n)

    E_old = 0.0
    MIXING = 0.3 
    
    # === Configuration Average (CA) Coefficients ===
    # For the Configuration Average of any p^N shell, the effective 
    # angular exchange coefficient for the F^2 integral is exactly 2/25.
    c_fock = 2.0 / 25.0 
    
    println("Starting ROHF SCF loop for 3p^3 (Configuration Average)...")
    if verbose
        @printf("%-4s | %-14s | %-10s | %-8s\n", "Iter", "E_total (Ha)", "Delta E", "Time (s)")
        println("-"^78)
    end

    for iter in 1:150 
        t0 = time()

        # ==========================================
        # PHASE 1: The Expanded Closed Core
        # ==========================================
        core_orbs = [orbitals[1], orbitals[2], orbitals[3], orbitals[4]]
        
        build_total_J_matrix!(ws, core_orbs)
        J_core .= ws.J
        
        assemble_K_matrix!(ws, ws.K_mats[0], 0, core_orbs)
        K_core_s .= ws.K_mats[0]
        
        assemble_K_matrix!(ws, ws.K_mats[1], 1, core_orbs)
        K_core_p .= ws.K_mats[1]

        # ==========================================
        # PHASE 2: The Open 3p Shell Contributions
        # ==========================================
        
        # A. Inter-shell: What the closed shells feel from the 3p shell (Occupancy = 3.0)
        build_specific_J_matrix!(ws, J_3p_spherical, orbitals[5], 3.0)
        
        assemble_K_matrix!(ws, ws.K_mats[0], 0, [orbitals[5]]) 
        K_3p_on_s .= ws.K_mats[0]
        
        assemble_K_matrix!(ws, ws.K_mats[1], 1, [orbitals[5]]) 
        K_3p_on_p .= ws.K_mats[1]

        # B. Intra-shell: What a 3p electron feels from its own shell (N-1 = 2.0)
        build_specific_J_matrix!(ws, J_3p_intra, orbitals[5], 2.0)
        
        # C. The Exact Exchange Quadrupole 
        # Inject the exact scalar for the Configuration Average
        build_specific_K_matrix!(ws, K_3p_intra, orbitals[5], 2, c_fock)

        # ==========================================
        # PHASE 3: Assemble Shell-Dependent Fock Matrices
        # ==========================================
        F_s = ws.F_mats[0]
        F_s .= H_core_s .+ J_core .+ J_3p_spherical .- K_core_s .- K_3p_on_s
        F_2p .= H_core_p .+ J_core .+ J_3p_spherical .- K_core_p .- K_3p_on_p
        F_3p .= H_core_p .+ J_core .+ J_3p_intra .- K_core_p .- K_3p_intra
        
        # ==========================================
        # PHASE 4: Diagonalization & Gram-Schmidt Projection
        # ==========================================
        evals_fs, evecs_fs = eigen(Symmetric(F_s[active_s, active_s]), ws.S[active_s, active_s])
        evals_f2p, evecs_f2p = eigen(Symmetric(F_2p[active_p, active_p]), ws.S[active_p, active_p])
        evals_f3p, evecs_f3p = eigen(Symmetric(F_3p[active_p, active_p]), ws.S[active_p, active_p])
        
        c_1s_new = zeros(Float64, n); c_1s_new[active_s] = evecs_fs[:, 1]
        c_2s_new = zeros(Float64, n); c_2s_new[active_s] = evecs_fs[:, 2]
        c_3s_new = zeros(Float64, n); c_3s_new[active_s] = evecs_fs[:, 3]
        
        c_2p_new = zeros(Float64, n); c_2p_new[active_p] = evecs_f2p[:, 1]
        c_3p_new = zeros(Float64, n); c_3p_new[active_p] = evecs_f3p[:, 2]
        
        # Gram-Schmidt Orthogonalization for the p-block
        overlap_2p_3p = dot(c_2p_new, ws.S * c_3p_new)
        c_3p_new .-= overlap_2p_3p .* c_2p_new
        
        c_1s_new ./= sqrt(dot(c_1s_new, ws.S * c_1s_new))
        c_2s_new ./= sqrt(dot(c_2s_new, ws.S * c_2s_new))
        c_3s_new ./= sqrt(dot(c_3s_new, ws.S * c_3s_new))
        c_2p_new ./= sqrt(dot(c_2p_new, ws.S * c_2p_new))
        c_3p_new ./= sqrt(dot(c_3p_new, ws.S * c_3p_new))
        
        orbitals[1].coeffs = MIXING * orbitals[1].coeffs + (1 - MIXING) * c_1s_new
        orbitals[2].coeffs = MIXING * orbitals[2].coeffs + (1 - MIXING) * c_2s_new
        orbitals[3].coeffs = MIXING * orbitals[3].coeffs + (1 - MIXING) * c_2p_new
        orbitals[4].coeffs = MIXING * orbitals[4].coeffs + (1 - MIXING) * c_3s_new
        orbitals[5].coeffs = MIXING * orbitals[5].coeffs + (1 - MIXING) * c_3p_new
        
        for orb in orbitals
            orb.coeffs ./= sqrt(dot(orb.coeffs, ws.S * orb.coeffs))
        end

        orbitals[1].energy = evals_fs[1]
        orbitals[2].energy = evals_fs[2]
        orbitals[3].energy = evals_f2p[1]
        orbitals[4].energy = evals_fs[3]
        orbitals[5].energy = evals_f3p[2] 
 
        # ==========================================
        # PHASE 5: Compute Total Energy
        # ==========================================
        E_total = 0.0

        for orb in orbitals
            if orb.occ > 0.0
                h_core = (orb.l == 0) ? H_core_s : H_core_p
                h_ii = dot(orb.coeffs, h_core * orb.coeffs)
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
            @printf("Final ROHF Energy (Config Avg): %.6f Ha\n", E_total)

            # for SO coupling later
            dense_grid = exp.(range(log(1e-8), log(R_max), length=10000))
            V_eff = compute_effective_central_potential(ws, orbitals, dense_grid, Z)
            P_3p = evaluate_orbital(ws.basis, orbitals[5].coeffs, dense_grid)
            
            filename = "phosphorus_rohf_CA_results_R$(R_max).jld2"
            jldsave(filename;
                    orbitals = orbitals,
                    E_total = E_total,
                    R_max = R_max,
                    V_nuclear = ws.V,
                    R_grid = dense_grid,
                    V_eff = V_eff,
                    P_3p = P_3p,
                    num_splines = n,
                    active_s = active_s,
                    active_p = active_p,
                    term = "ConfigAvg"
                    )
            println("Saved Configuration Average ROHF data to $filename")
            println("===== END =====")
            break
        end
        
        E_old = E_total
    end
end
function compute_effective_central_potential(ws::SolverWorkspace, orbitals::Vector{Orbital}, dense_grid::Vector{Float64}, Z::Float64)
    n_splines = ws.basis.num_splines
    total_Y0_coeffs = zeros(Float64, n_splines)
    y_temp = zeros(Float64, n_splines)
    
    # 1. Superposition of the Hartree Screening Functions
    for orb in orbitals
        if orb.occ > 0.0
            # Solve the Poisson equation for the individual orbital's density
            solve_poisson_J!(ws, y_temp, orb)
            
            # Linearly scale by fractional occupancy and accumulate
            total_Y0_coeffs .+= orb.occ .* y_temp
        end
    end
    
    # 2. Evaluate the continuous spline over the dense physical grid
    Y0_grid = evaluate_orbital(ws.basis, total_Y0_coeffs, dense_grid)
    
    # 3. Construct V_eff(r) = (-Z + Y0(r)) / r
    n_points = length(dense_grid)
    V_eff = zeros(Float64, n_points)
    
    for i in 1:n_points
        r = dense_grid[i]
        # Prevent division by zero precisely at the origin
        if r > 1e-12
            V_eff[i] = (-Z + Y0_grid[i]) / r
        else
            # Deep in the core, the potential is heavily dominated by the bare nucleus
            V_eff[i] = -Z / 1e-12 
        end
    end
    
    return V_eff
end

solve_phosphorus_rohf_CA(30.0)
