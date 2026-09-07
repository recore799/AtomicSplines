using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))

using AtomicSplines
using LinearAlgebra
using Printf
using JLD2

isdefined(Main, :DIISState) || include(joinpath(@__DIR__, "rohf_diis.jl"))

"""
    solve_tin_rohf(R_max; estado, use_diis, tol, max_iter, diis_thresh, save)

`estado = nothing` conserva el prompt interactivo. Con `use_diis` el level shift deja
de ser permanente: se mantiene mientras el residual del conmutador siga por encima de
`diis_thresh` y se apaga en cuanto se enciende la extrapolacion de Pulay.
"""
function solve_tin_rohf(R_max; verbose::Bool=true, estado=nothing,
                        use_diis::Bool=false, tol::Float64=1e-9,
                        max_iter::Int=1000, diis_thresh::Float64=1e-2,
                        save::Bool=true,
                        coeff_k2_override::Union{Nothing,Float64}=nothing)
    if estado === nothing
        print("¿A qué estado desea optimizar? (av / 3P): ")
        estado = strip(readline())
    end
    if estado ∉ ["av", "3P"]
        println("Estado no reconocido. Usando 'av' por defecto.")
        estado = "av"
    end
    coeff_k2 = (estado == "av") ? (2.0 / 25.0) : (5.0 / 25.0)
    # Solo para verificacion: permite forzar el coeficiente de intercambio intra-capa
    # sin tocar la fisica por defecto. El valor por defecto (nothing) no cambia nada.
    coeff_k2_override === nothing || (coeff_k2 = coeff_k2_override)

    println("=== Tin ROHF Optimización: $estado (Z=50) ===")
    
    N_elems = 500
    Z = 50.0

    ws = cached_init_scf_workspace(R_max, N_elems, Val(8), Z; γ=4.0, calc_R_matrices=true)
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
        Orbital(3, 2, 10.0), # [6] 3d (Closed)
        Orbital(4, 0, 2.0),  # [7] 4s (Closed)
        Orbital(4, 1, 6.0),  # [8] 4p (Closed)
        Orbital(4, 2, 10.0), # [9] 4d (Closed)
        Orbital(5, 0, 2.0),  # [10] 5s (Closed)
        Orbital(5, 1, 2.0)   # [11] 5p (Open Valence)
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
    orbitals[10].coeffs = zeros(Float64, n); orbitals[10].coeffs[active_s] = evecs_s[:, 5]
    
    evals_p, evecs_p = eigen(Symmetric(H_core_p[active_p, active_p]), ws.S[active_p, active_p])
    orbitals[3].coeffs = zeros(Float64, n); orbitals[3].coeffs[active_p] = evecs_p[:, 1]
    orbitals[5].coeffs = zeros(Float64, n); orbitals[5].coeffs[active_p] = evecs_p[:, 2]
    orbitals[8].coeffs = zeros(Float64, n); orbitals[8].coeffs[active_p] = evecs_p[:, 3]
    orbitals[11].coeffs = zeros(Float64, n); orbitals[11].coeffs[active_p] = evecs_p[:, 4]

    evals_d, evecs_d = eigen(Symmetric(H_core_d[active_d, active_d]), ws.S[active_d, active_d])
    orbitals[6].coeffs = zeros(Float64, n); orbitals[6].coeffs[active_d] = evecs_d[:, 1]
    orbitals[9].coeffs = zeros(Float64, n); orbitals[9].coeffs[active_d] = evecs_d[:, 2]
    
    for orb in orbitals
        orb.coeffs ./= sqrt(dot(orb.coeffs, ws.S * orb.coeffs))
    end

    # Pre-allocate temporary matrices
    J_core = zeros(Float64, n, n); K_core_s = zeros(Float64, n, n)
    K_core_p = zeros(Float64, n, n); K_core_d = zeros(Float64, n, n)
    J_5p_spherical = zeros(Float64, n, n); K_5p_on_s = zeros(Float64, n, n)
    K_5p_on_p = zeros(Float64, n, n); K_5p_on_d = zeros(Float64, n, n)
    J_5p_intra = zeros(Float64, n, n); K_5p_intra = zeros(Float64, n, n)
    
    F_core_p = zeros(Float64, n, n); F_core_d = zeros(Float64, n, n); F_5p = zeros(Float64, n, n)

    E_old = 0.0
    MIXING = 0.0
    diis = DIISState(; thresh = diis_thresh)
    diis_on = false
    # (iter, E_total, |dE|, residual proyectado, conmutador crudo) por iteracion:
    # es la curva que va a la figura de convergencia de resultados.tex.
    trace = NTuple{5,Float64}[]
    
    println("Comenzando ciclo SCF ROHF...")
    if verbose
        @printf("%-4s | %-14s | %-10s | %-8s | %-10s | %s\n",
                "Iter", "E_total (Ha)", "Delta E", "Time (s)", "|FDS-SDF|", "modo")
        println("-"^78)
    end

    for iter in 1:max_iter
        t0 = time()

        # --- Phase 1: Closed Core Matrix Assembly ---
        core_orbs = orbitals[1:10]
        build_total_J_matrix!(ws, core_orbs); J_core .= ws.J
        assemble_K_matrix!(ws, ws.K_mats[0], 0, core_orbs); K_core_s .= ws.K_mats[0]
        assemble_K_matrix!(ws, ws.K_mats[1], 1, core_orbs); K_core_p .= ws.K_mats[1]
        assemble_K_matrix!(ws, ws.K_mats[2], 2, core_orbs); K_core_d .= ws.K_mats[2]

        # --- Phase 2: Open 5p Shell Contributions ---
        build_specific_J_matrix!(ws, J_5p_spherical, orbitals[11], 2.0)
        assemble_K_matrix!(ws, ws.K_mats[0], 0, [orbitals[11]]); K_5p_on_s .= ws.K_mats[0]
        assemble_K_matrix!(ws, ws.K_mats[1], 1, [orbitals[11]]); K_5p_on_p .= ws.K_mats[1]
        assemble_K_matrix!(ws, ws.K_mats[2], 2, [orbitals[11]]); K_5p_on_d .= ws.K_mats[2]

        build_specific_J_matrix!(ws, J_5p_intra, orbitals[11], 1.0)
        build_specific_K_matrix!(ws, K_5p_intra, orbitals[11], 2, coeff_k2)

        # --- Phase 3: Raw Fock Matrices ---
        F_s = ws.F_mats[0]
        F_s .= H_core_s .+ J_core .+ J_5p_spherical .- K_core_s .- K_5p_on_s
        F_core_p .= H_core_p .+ J_core .+ J_5p_spherical .- K_core_p .- K_5p_on_p
        F_core_d .= H_core_d .+ J_core .+ J_5p_spherical .- K_core_d .- K_5p_on_d
        F_5p .= H_core_p .+ J_core .+ J_5p_intra .- K_core_p .- K_5p_intra

        # ==========================================
        # PHASE 3b: C-DIIS por bloques (ver rohf_diis.jl)
        # ==========================================
        # Cuatro canales, cada uno con su pareja (F, D): las s cerradas, las p de core
        # (2p, 3p, 4p), las d de core (3d, 4d) y el 5p abierto con su ocupacion
        # fraccionaria. F_core_p y F_5p comparten espacio activo pero son operadores
        # distintos, asi que sus densidades no se pueden intercambiar. El historial
        # guarda las Fock CRUDAS: el level shift se aplica despues, y solo mientras
        # DIIS siga apagado.
        F_use = Matrix{Float64}[F_s, F_core_p, F_core_d, F_5p]
        if use_diis
            D_s  = block_density(orbitals, (1, 2, 4, 7, 10), active_s)
            D_cp = block_density(orbitals, (3, 5, 8), active_p)
            D_cd = block_density(orbitals, (6, 9), active_d)
            D_5p = block_density(orbitals, (11,), active_p)
            # El espacio ocupado l=1 incluye al 5p abierto: es contra el que se
            # ortogonaliza el core y de donde salen los multiplicadores de Lagrange.
            C_s  = occ_block(orbitals, (1, 2, 4, 7, 10), active_s)
            C_p  = occ_block(orbitals, (3, 5, 8, 11), active_p)
            C_d  = occ_block(orbitals, (6, 9), active_d)
            F_use, diis_on = diis_update!(diis, iter,
                                          Matrix{Float64}[F_s, F_core_p, F_core_d, F_5p],
                                          Matrix{Float64}[D_s, D_cp, D_cd, D_5p],
                                          Matrix{Float64}[C_s, C_p, C_d, C_p],
                                          UnitRange{Int}[active_s, active_p, active_d, active_p],
                                          ws.S)
        end
        F_s_use, F_cp_use, F_cd_use, F_5p_use = F_use[1], F_use[2], F_use[3], F_use[4]

        # ==========================================
        # PHASE 4: Scheduled Dynamic Level Shifting
        # ==========================================
        # Un desplazamiento constante de 3 Ha que nunca se apaga amortigua la
        # convergencia de forma permanente. Con C-DIIS el shift es solo el regimen de
        # arranque y se suelta al entrar en la cuenca de atraccion.
        level_shift = diis_on ? 0.0 : 3.0

        if level_shift > 0.0
            # Construct purely spatial geometric projectors
            P_s_spatial = zeros(Float64, n, n)
            for i in [1, 2, 4, 7, 10] # 1s, 2s, 3s, 4s, 5s
                P_s_spatial .+= orbitals[i].coeffs * orbitals[i].coeffs'
            end
            
            P_p_spatial = zeros(Float64, n, n)
            for i in [3, 5, 8, 11] # 2p, 3p, 4p, 5p
                P_p_spatial .+= orbitals[i].coeffs * orbitals[i].coeffs'
            end
            
            P_d_spatial = zeros(Float64, n, n)
            for i in [6, 9] # 3d, 4d
                P_d_spatial .+= orbitals[i].coeffs * orbitals[i].coeffs'
            end

            # Compute covariant projection matrices
            S_Ps_S = ws.S * P_s_spatial * ws.S
            S_Pp_S = ws.S * P_p_spatial * ws.S
            S_Pd_S = ws.S * P_d_spatial * ws.S

            # El shift se aplica sobre COPIAS: las Fock crudas se necesitan intactas
            # para las energias orbitales de la fase 6.
            F_s_use  = copy(F_s_use);  F_cp_use = copy(F_cp_use)
            F_cd_use = copy(F_cd_use); F_5p_use = copy(F_5p_use)

            # Apply shift to the virtual space
            F_s_use  .+= level_shift .* (ws.S .- S_Ps_S)
            F_cp_use .+= level_shift .* (ws.S .- S_Pp_S)
            F_cd_use .+= level_shift .* (ws.S .- S_Pd_S)
            F_5p_use .+= level_shift .* (ws.S .- S_Pp_S)
            
            if verbose && iter <= 21
                @printf("  [Level Shift Active: b = %.3f Ha]\n", level_shift)
            end
        end

        # --- Phase 5: Diagonalization & Gram-Schmidt ---
        evals_fs, evecs_fs = eigen(Symmetric(F_s_use[active_s, active_s]), ws.S[active_s, active_s])
        evals_fcp, evecs_fcp = eigen(Symmetric(F_cp_use[active_p, active_p]), ws.S[active_p, active_p])
        evals_fcd, evecs_fcd = eigen(Symmetric(F_cd_use[active_d, active_d]), ws.S[active_d, active_d])
        evals_f5p, evecs_f5p = eigen(Symmetric(F_5p_use[active_p, active_p]), ws.S[active_p, active_p])
        
        c_1s_new = zeros(Float64, n); c_1s_new[active_s] = evecs_fs[:, 1]
        c_2s_new = zeros(Float64, n); c_2s_new[active_s] = evecs_fs[:, 2]
        c_3s_new = zeros(Float64, n); c_3s_new[active_s] = evecs_fs[:, 3]
        c_4s_new = zeros(Float64, n); c_4s_new[active_s] = evecs_fs[:, 4]
        c_5s_new = zeros(Float64, n); c_5s_new[active_s] = evecs_fs[:, 5]
        
        c_2p_new = zeros(Float64, n); c_2p_new[active_p] = evecs_fcp[:, 1]
        c_3p_new = zeros(Float64, n); c_3p_new[active_p] = evecs_fcp[:, 2]
        c_4p_new = zeros(Float64, n); c_4p_new[active_p] = evecs_fcp[:, 3]
        c_5p_new = zeros(Float64, n); c_5p_new[active_p] = evecs_f5p[:, 4]
        
        c_3d_new = zeros(Float64, n); c_3d_new[active_d] = evecs_fcd[:, 1]
        c_4d_new = zeros(Float64, n); c_4d_new[active_d] = evecs_fcd[:, 2]
        
        # Cascaded Gram-Schmidt Orthogonalization
        c_2s_new .-= dot(c_1s_new, ws.S * c_2s_new) .* c_1s_new
        c_3s_new .-= dot(c_1s_new, ws.S * c_3s_new) .* c_1s_new
        c_3s_new .-= dot(c_2s_new, ws.S * c_3s_new) .* c_2s_new
        c_4s_new .-= dot(c_1s_new, ws.S * c_4s_new) .* c_1s_new
        c_4s_new .-= dot(c_2s_new, ws.S * c_4s_new) .* c_2s_new
        c_4s_new .-= dot(c_3s_new, ws.S * c_4s_new) .* c_3s_new
        c_5s_new .-= dot(c_1s_new, ws.S * c_5s_new) .* c_1s_new
        c_5s_new .-= dot(c_2s_new, ws.S * c_5s_new) .* c_2s_new
        c_5s_new .-= dot(c_3s_new, ws.S * c_5s_new) .* c_3s_new
        c_5s_new .-= dot(c_4s_new, ws.S * c_5s_new) .* c_4s_new
        
        c_3p_new .-= dot(c_2p_new, ws.S * c_3p_new) .* c_2p_new
        c_4p_new .-= dot(c_2p_new, ws.S * c_4p_new) .* c_2p_new
        c_4p_new .-= dot(c_3p_new, ws.S * c_4p_new) .* c_3p_new
        c_5p_new .-= dot(c_2p_new, ws.S * c_5p_new) .* c_2p_new
        c_5p_new .-= dot(c_3p_new, ws.S * c_5p_new) .* c_3p_new
        c_5p_new .-= dot(c_4p_new, ws.S * c_5p_new) .* c_4p_new
        
        c_4d_new .-= dot(c_3d_new, ws.S * c_4d_new) .* c_3d_new
        
        c_1s_new ./= sqrt(dot(c_1s_new, ws.S * c_1s_new)); c_2s_new ./= sqrt(dot(c_2s_new, ws.S * c_2s_new))
        c_3s_new ./= sqrt(dot(c_3s_new, ws.S * c_3s_new)); c_4s_new ./= sqrt(dot(c_4s_new, ws.S * c_4s_new))
        c_5s_new ./= sqrt(dot(c_5s_new, ws.S * c_5s_new))
        
        c_2p_new ./= sqrt(dot(c_2p_new, ws.S * c_2p_new)); c_3p_new ./= sqrt(dot(c_3p_new, ws.S * c_3p_new))
        c_4p_new ./= sqrt(dot(c_4p_new, ws.S * c_4p_new)); c_5p_new ./= sqrt(dot(c_5p_new, ws.S * c_5p_new))
        
        c_3d_new ./= sqrt(dot(c_3d_new, ws.S * c_3d_new)); c_4d_new ./= sqrt(dot(c_4d_new, ws.S * c_4d_new))
        
        # Linear Damping
        orbitals[1].coeffs = MIXING * orbitals[1].coeffs + (1 - MIXING) * c_1s_new
        orbitals[2].coeffs = MIXING * orbitals[2].coeffs + (1 - MIXING) * c_2s_new
        orbitals[3].coeffs = MIXING * orbitals[3].coeffs + (1 - MIXING) * c_2p_new
        orbitals[4].coeffs = MIXING * orbitals[4].coeffs + (1 - MIXING) * c_3s_new
        orbitals[5].coeffs = MIXING * orbitals[5].coeffs + (1 - MIXING) * c_3p_new
        orbitals[6].coeffs = MIXING * orbitals[6].coeffs + (1 - MIXING) * c_3d_new
        orbitals[7].coeffs = MIXING * orbitals[7].coeffs + (1 - MIXING) * c_4s_new
        orbitals[8].coeffs = MIXING * orbitals[8].coeffs + (1 - MIXING) * c_4p_new
        orbitals[9].coeffs = MIXING * orbitals[9].coeffs + (1 - MIXING) * c_4d_new
        orbitals[10].coeffs = MIXING * orbitals[10].coeffs + (1 - MIXING) * c_5s_new
        orbitals[11].coeffs = MIXING * orbitals[11].coeffs + (1 - MIXING) * c_5p_new
        
        for orb in orbitals; orb.coeffs ./= sqrt(dot(orb.coeffs, ws.S * orb.coeffs)); end

        orbitals[1].energy = evals_fs[1]; orbitals[2].energy = evals_fs[2]; orbitals[3].energy = evals_fcp[1]
        orbitals[4].energy = evals_fs[3]; orbitals[5].energy = evals_fcp[2]; orbitals[6].energy = evals_fcd[1]
        orbitals[7].energy = evals_fs[4]; orbitals[8].energy = evals_fcp[3]; orbitals[9].energy = evals_fcd[2]
        orbitals[10].energy = evals_fs[5]; orbitals[11].energy = evals_f5p[4]
 
        # --- Phase 6: Compute Total Energy ---
        E_total = 0.0
        for orb in orbitals
            if orb.occ > 0.0
                h_core = orb.l == 0 ? H_core_s : (orb.l == 1 ? H_core_p : H_core_d)
                E_total += (orb.occ / 2.0) * (dot(orb.coeffs, h_core * orb.coeffs) + orb.energy)
            end
        end
        
        delta = abs(E_total - E_old)
        push!(trace, (Float64(iter), E_total, delta, diis.resid, diis.resid_raw))
        elapsed = time() - t0
        
        if verbose
            @printf("%-4d | %14.8f | %10.2e | %8.4f | %10.2e | %s\n",
                    iter, E_total, delta, elapsed, diis.resid, diis.frozen ? "libre" : (diis_on ? "DIIS" : "shift"))
        end
        
        # Pure Energy Stabilization Metric
        # No se acepta la convergencia mientras DIIS este extrapolando con el residual
        # ESTANCADO: ahi |dE| < tol puede dispararse sobre un ciclo limite en vez de
        # sobre un punto fijo. Con el residual todavia bajando (stall == 0) el minimo
        # es autentico y se acepta igual que sin DIIS.
        if delta < tol && (!diis.on || diis.frozen || diis.stall == 0)
            if verbose 
                println("-"^78)
                println("Converged purely on energy metric in $iter iterations.")
            end
            @printf("Energía final HF (%s): %.6f Ha\n", estado, E_total)

            if !save
                return (E_total = E_total, iters = iter, converged = true,
                        switch_iter = diis.switch_iter, freeze_iter = diis.freeze_iter,
                        restarts = diis.restarts,
                        resid = diis.resid, trace = trace, orbitals = orbitals)
            end

            dense_grid = exp.(range(log(1e-8), log(R_max), length=10000))
            V_eff = compute_effective_central_potential(ws, orbitals, dense_grid, Z)
            P_5p = evaluate_orbital(ws.basis, orbitals[11].coeffs, dense_grid)

            filename = "tin_rohf_results_$(estado)_R$(R_max).jld2"
            
            jldsave(filename;
                orbitals = orbitals, E_total = E_total, R_max = R_max,
                R_grid = dense_grid, 
                V_nuclear = ws.V, num_splines = n, active_s = active_s,
                active_p = active_p, active_d = active_d,
                V_eff = V_eff, P_5p = P_5p
            )
            println("Saved Term-Dependent ROHF data and SO coupling prerequisites to $filename")
            println("===== END =====")
            return (E_total = E_total, iters = iter, converged = true,
                    switch_iter = diis.switch_iter, freeze_iter = diis.freeze_iter,
                        restarts = diis.restarts,
                    resid = diis.resid, trace = trace, orbitals = orbitals)
        end
        
        E_old = E_total
    end
    @warn "Estanio: no convergio en $max_iter iteraciones (tol = $tol)."
    return (E_total = E_old, iters = max_iter, converged = false,
            switch_iter = diis.switch_iter, freeze_iter = diis.freeze_iter,
                        restarts = diis.restarts,
            resid = diis.resid, trace = trace, orbitals = orbitals)
end

if abspath(PROGRAM_FILE) == @__FILE__
    # Execute the default 3P state directly if run from script
    solve_tin_rohf(30.0)
end
