using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines

using LinearAlgebra
using Printf
using JLD2

function solve_lithium_rohf(R_max; verbose::Bool=true)
    println("=== Lithium ROHF Term-Dependent 2S with C-DIIS (Z=3) ===")
    
    N_elems = 200
    Z = 3.0
    basis = generate_basis(R_max, N_elems, Val(7), γ=3.0)
    ws = init_scf_workspace(basis, Z)

    n = basis.num_splines
    active_s = 2:(n-1)  
    
    orbitals = [
        Orbital(1, 0, 2.0), # 1s
        Orbital(2, 0, 1.0)  # 2s
    ]
    
    H_core_s = ws.T + ws.V
    
    # --- Initial Guess ---
    evals_s, evecs_s = eigen(Symmetric(H_core_s[active_s, active_s]), ws.S[active_s, active_s])
    orbitals[1].coeffs = zeros(Float64, n); orbitals[1].coeffs[active_s] = evecs_s[:, 1]
    orbitals[2].coeffs = zeros(Float64, n); orbitals[2].coeffs[active_s] = evecs_s[:, 2]
    
    for orb in orbitals
        orb.coeffs ./= sqrt(dot(orb.coeffs, ws.S * orb.coeffs))
    end

    J_1s = zeros(Float64, n, n)
    K_1s = zeros(Float64, n, n)
    J_2s = zeros(Float64, n, n)
    K_2s = zeros(Float64, n, n)
    
    F_1s = zeros(Float64, n, n)
    F_2s = zeros(Float64, n, n)
    F_U  = zeros(Float64, n, n)

    # --- DIIS Subspace Initialization ---
    DIIS_MAX_VECS = 15
    F1s_hist = Matrix{Float64}[]
    F2s_hist = Matrix{Float64}[]
    err_hist = Vector{Float64}[]

    E_old = 0.0
    
    println("Comenzando ciclo SCF ROHF con C-DIIS para Litio...")
    if verbose
        @printf("%-4s | %-14s | %-10s | %-10s | %-8s\n", "Iter", "E_total (Ha)", "Delta E", "Max Error", "Time (s)")
        println("-"^78)
    end

    for iter in 1:500
        t0 = time()

        # 1. Build Shell Potentials
        build_specific_J_matrix!(ws, J_1s, orbitals[1], 2.0)
        assemble_K_matrix!(ws, ws.K_mats[0], 0, [orbitals[1]])
        K_1s .= ws.K_mats[0]
        
        build_specific_J_matrix!(ws, J_2s, orbitals[2], 1.0)
        assemble_K_matrix!(ws, ws.K_mats[0], 0, [orbitals[2]])
        K_2s .= ws.K_mats[0]

        F_1s .= H_core_s .+ J_1s .+ J_2s .- K_1s .- K_2s
        F_2s .= H_core_s .+ J_1s .- K_1s
        
        # 2. Compute Density Matrices
        D_1s = orbitals[1].coeffs * orbitals[1].coeffs'
        D_2s = orbitals[2].coeffs * orbitals[2].coeffs'
        D_tot = D_1s .+ D_2s
        
        # 3. Construct the Raw Unified Matrix
        Delta_F = F_1s .- F_2s
        term1 = ws.S * D_1s * Delta_F
        term2 = Delta_F * D_1s * ws.S
        term3 = ws.S * D_2s * Delta_F
        term4 = Delta_F * D_2s * ws.S
        
        F_U .= F_1s .- term1 .- term2 .+ term3 .+ term4
        F_U .= Symmetric(0.5 .* (F_U .+ F_U'))
        
        # 4. Compute the C-DIIS Error Commutator Vector
        err_mat = F_U * D_tot * ws.S .- ws.S * D_tot * F_U
        err_vec = vec(err_mat)
        max_err = maximum(abs.(err_vec))
        
        # Update Subspace Histories
        push!(F1s_hist, copy(F_1s))
        push!(F2s_hist, copy(F_2s))
        push!(err_hist, copy(err_vec))
        
        if length(err_hist) > DIIS_MAX_VECS
            popfirst!(F1s_hist)
            popfirst!(F2s_hist)
            popfirst!(err_hist)
        end
        
        # 5. Execute Pulay Extrapolation (Only if subspace has accumulated)
        num_vecs = length(err_hist)
        if num_vecs > 1
            # Build the symmetric Pulay B matrix
            B = zeros(Float64, num_vecs + 1, num_vecs + 1)
            for i in 1:num_vecs
                for j in i:num_vecs
                    val = dot(err_hist[i], err_hist[j])
                    B[i, j] = val
                    B[j, i] = val
                end
                B[i, num_vecs + 1] = -1.0
                B[num_vecs + 1, i] = -1.0
            end
            B[num_vecs + 1, num_vecs + 1] = 0.0
            
            # Add microscopic dampening to prevent singular exceptions
            for i in 1:num_vecs
                B[i, i] += 1e-12
            end
            
            # Build the RHS vector [-0, -0, ... , -1]
            rhs = zeros(Float64, num_vecs + 1)
            rhs[num_vecs + 1] = -1.0
            
            # Solve for the optimal extrapolation weights
            weights = B \ rhs
            
            # Extrapolate the Physical Potentials
            F_1s_ext = zeros(Float64, n, n)
            F_2s_ext = zeros(Float64, n, n)
            for i in 1:num_vecs
                F_1s_ext .+= weights[i] .* F1s_hist[i]
                F_2s_ext .+= weights[i] .* F2s_hist[i]
            end
            
            # Re-couple the Extrapolated Unified Matrix
            Delta_F_ext = F_1s_ext .- F_2s_ext
            term1_ext = ws.S * D_1s * Delta_F_ext
            term2_ext = Delta_F_ext * D_1s * ws.S
            term3_ext = ws.S * D_2s * Delta_F_ext
            term4_ext = Delta_F_ext * D_2s * ws.S
            
            F_U .= F_1s_ext .- term1_ext .- term2_ext .+ term3_ext .+ term4_ext
            F_U .= Symmetric(0.5 .* (F_U .+ F_U'))
        end
        
        # 6. Diagonalize the Stabilized Operator
        evals_fu, evecs_fu = eigen(Symmetric(F_U[active_s, active_s]), ws.S[active_s, active_s])
        
        c_1s_new = zeros(Float64, n); c_1s_new[active_s] = evecs_fu[:, 1]
        c_2s_new = zeros(Float64, n); c_2s_new[active_s] = evecs_fu[:, 2]
        
        c_1s_new ./= sqrt(dot(c_1s_new, ws.S * c_1s_new))
        c_2s_new ./= sqrt(dot(c_2s_new, ws.S * c_2s_new))
        
        # Immediate overwrite (DIIS completely replaces the need for Linear Mixing)
        orbitals[1].coeffs .= c_1s_new
        orbitals[2].coeffs .= c_2s_new
        
        orbitals[1].energy = dot(orbitals[1].coeffs, F_1s * orbitals[1].coeffs)
        orbitals[2].energy = dot(orbitals[2].coeffs, F_2s * orbitals[2].coeffs)
 
        # 7. Compute Energy
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
            @printf("%-4d | %14.8f | %10.2e | %10.2e | %8.4f\n", iter, E_total, delta, max_err, elapsed)
        end
        
        # Convergence is now strictly gated by the physical density commutator reaching zero
        if max_err < 1e-8 && delta < 1e-9 && iter > 2
            if verbose 
                println("-"^78)
                println("C-DIIS Converged perfectly in $iter iterations.")
            end
            @printf("Energía final HF-t (^2S): %.8f Ha\n", E_total)

            filename = "lithium_rohf_cdiis_R$(R_max).jld2"
            jldsave(filename;
                orbitals = orbitals,
                E_total = E_total,
                R_max = R_max,
                R_grid = ws.R,
                V_nuclear = ws.V,
                num_splines = n,
                active_s = active_s
            )
            println("Saved DIIS optimized Unified ROHF data to $filename")
            println("===== END =====")
            break
        end
        
        E_old = E_total
    end
end

solve_lithium_rohf(30.0)
