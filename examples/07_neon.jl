using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using LinearAlgebra
using Printf

function solve_neon(R_max)
    println("=== Neon (Z=10) ===")
    
    N_elems = 100
    Z = 10.0
    basis = generate_basis(R_max, N_elems, Val(7), γ=2.5)
    ws = init_scf_workspace(basis, Z)

    # Clear cached Poisson factors from previous runs if workspace is reused NEEDS FIXING
    empty!(ws.poisson_factors)

    n = basis.num_splines

    # Enforce origin boundary conditions: P(r) ~ r^(l+1)
    # s-orbitals (l=0): Drop B1 to enforce P(0) = 0
    active_s = 2:(n-1) 
    # p-orbitals (l=1): Drop B1 and B2 to enforce P(0) = 0 and P'(0) = 0
    active_p = 3:(n-1) 
    
    # Initialize Orbitals (n, l, occupancy)
    orbitals = [
        Orbital(1, 0, 2.0), # 1s
        Orbital(2, 0, 2.0), # 2s
        Orbital(2, 1, 6.0)  # 2p
    ]
    
    # --- Core Hamiltonians ---
    # s-block (l=0): Kinetic + Electron-Nuclear
    H_core_s = ws.T + ws.V
    
    # p-block (l=1): Kinetic + Electron-Nuclear + Centrifugal (1/r^2 = ws.V2)
    H_core_p = ws.T + ws.V + ws.V2 
    
    # --- Initial Guess ---
    # Diagonalize the core Hamiltonian to get starting coefficients
    evals_s, evecs_s = eigen(Matrix(H_core_s[active_s, active_s]), Matrix(ws.S[active_s, active_s]))
    orbitals[1].coeffs = zeros(Float64, n); orbitals[1].coeffs[active_s] = evecs_s[:, 1]
    orbitals[2].coeffs = zeros(Float64, n); orbitals[2].coeffs[active_s] = evecs_s[:, 2]
    
    evals_p, evecs_p = eigen(Matrix(H_core_p[active_p, active_p]), Matrix(ws.S[active_p, active_p]))
    orbitals[3].coeffs = zeros(Float64, n); orbitals[3].coeffs[active_p] = evecs_p[:, 1]
    
    # Normalize initial guesses
    for orb in orbitals
        orb.coeffs ./= sqrt(dot(orb.coeffs, ws.S * orb.coeffs))
    end

    # --- SCF Loop ---
    E_old = 0.0
    MIXING = 0.0 # Dampening factor to aid convergence
    
    println("Comenzando ciclo SCF...")
    
    for iter in 1:60
        t0 = time()
        
        # --- Build Potentials ---
        y_total = zeros(Float64, n)
        for orb in orbitals
            y_orb = solve_poisson_J(ws, orb)
            y_total .+= orb.occ .* y_orb
        end
        
        # Coulomb matrix is isotropic (k=0) for closed shells
        J_mat = assemble_J_matrix(ws, y_total)
        
        # Exchange matrices are l-dependent
        K_s = assemble_K_matrix(ws, 0, orbitals)
        K_p = assemble_K_matrix(ws, 1, orbitals)
        
        # --- Build Fock Matrices ---
        F_s = H_core_s + J_mat - K_s
        F_p = H_core_p + J_mat - K_p
        
        # --- Diagonalize Independent Blocks ---
        evals_fs, evecs_fs = eigen(Matrix(F_s[active_s, active_s]), Matrix(ws.S[active_s, active_s]))
        evals_fp, evecs_fp = eigen(Matrix(F_p[active_p, active_p]), Matrix(ws.S[active_p, active_p]))
        
        # --- Update Orbitals (With Mixing) ---
        # 1s
        c_1s_new = zeros(Float64, n); c_1s_new[active_s] = evecs_fs[:, 1]
        c_1s_new ./= sqrt(dot(c_1s_new, ws.S * c_1s_new))
        if dot(orbitals[1].coeffs, ws.S * c_1s_new) < 0
            c_1s_new .= -c_1s_new # Enforce consistent phase
        end
        orbitals[1].coeffs = MIXING * orbitals[1].coeffs + (1 - MIXING) * c_1s_new
        orbitals[1].coeffs ./= sqrt(dot(orbitals[1].coeffs, ws.S * orbitals[1].coeffs))
        orbitals[1].energy = evals_fs[1]

        # 2s
        c_2s_new = zeros(Float64, n); c_2s_new[active_s] = evecs_fs[:, 2]
        c_2s_new ./= sqrt(dot(c_2s_new, ws.S * c_2s_new))
        if dot(orbitals[2].coeffs, ws.S * c_2s_new) < 0
            c_2s_new .= -c_2s_new
        end
        orbitals[2].coeffs = MIXING * orbitals[2].coeffs + (1 - MIXING) * c_2s_new
        orbitals[2].coeffs ./= sqrt(dot(orbitals[2].coeffs, ws.S * orbitals[2].coeffs))
        orbitals[2].energy = evals_fs[2]

        # 2p
        c_2p_new = zeros(Float64, n); c_2p_new[active_p] = evecs_fp[:, 1]
        c_2p_new ./= sqrt(dot(c_2p_new, ws.S * c_2p_new))
        if dot(orbitals[3].coeffs, ws.S * c_2p_new) < 0
            c_2p_new .= -c_2p_new
        end
        orbitals[3].coeffs = MIXING * orbitals[3].coeffs + (1 - MIXING) * c_2p_new
        orbitals[3].coeffs ./= sqrt(dot(orbitals[3].coeffs, ws.S * orbitals[3].coeffs))
        orbitals[3].energy = evals_fp[1]
 
        # --- Compute Total Energy ---
        # E = Σ_i (occ_i / 2) * (h_ii + ε_i)
        E_total = 0.0
        
        h_11 = dot(orbitals[1].coeffs, H_core_s * orbitals[1].coeffs)
        E_total += (orbitals[1].occ / 2.0) * (h_11 + orbitals[1].energy)
        
        h_22 = dot(orbitals[2].coeffs, H_core_s * orbitals[2].coeffs)
        E_total += (orbitals[2].occ / 2.0) * (h_22 + orbitals[2].energy)
        
        h_pp = dot(orbitals[3].coeffs, H_core_p * orbitals[3].coeffs)
        E_total += (orbitals[3].occ / 2.0) * (h_pp + orbitals[3].energy)
        
        delta = abs(E_total - E_old)
        
        if delta < 1e-9
            println("Radio de confinamiento: $R_max a.u.")
            @printf("Energía final: %.6f Ha\n", E_total)
            println("===== END =====")
            break
        end
        
        E_old = E_total
    end
end

# r_c_values = [1.0:0.35:1.35; 1.4:0.05:1.5; 1.75:0.05:1.85; 2.0:0.5:4.5]

# for r_c in r_c_values
#     solve_neon(r_c)
# end

solve_neon(10.0)
