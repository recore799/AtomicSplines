using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))

using AtomicSplines
using LinearAlgebra
using Printf
using JLD2

isdefined(Main, :DIISState) || include(joinpath(@__DIR__, "rohf_diis.jl"))

"""
    solve_silicon_rohf(R_max; estado, use_diis, tol, max_iter, diis_thresh, save)

`estado = nothing` conserva el prompt interactivo. `use_diis` enciende el C-DIIS por
bloques de `rohf_diis.jl` y `save = false` evita escribir el .jld2.
"""
function solve_silicon_rohf(R_max; verbose::Bool=true, estado=nothing,
                            use_diis::Bool=false, tol::Float64=1e-9,
                            max_iter::Int=150, diis_thresh::Float64=1e-2,
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
    # f_2(3P) = -5/25 por Slater-Condon (lo confirma el test T3 de np2_ci_full.jl: las
    # diferencias valen 0.24 y 0.60 F^2). El ensamblado de abajo RESTA K_intra, asi que
    # la convencion es coeff_k2 = -f_2 y el 3P exige +5/25, igual que en germanio y
    # estanio. Con el -5/25 anterior el 3P salia 0.094 Ha (C) / 0.065 Ha (Si) POR ENCIMA
    # del limite Hartree-Fock y por encima del promedio de configuracion, que es
    # imposible para el termino fundamental. Corregido reproduce a Froese Fischer y a
    # tab:resultados_energia_global del manuscrito a 8 cifras.
    # Los *_rohf_results_3P_R30.0.jld2 de C y Si ya se reemplazaron por los que produce
    # este signo; estan trackeados y registrados en benchmarks/RESULTS.toml.
    coeff_k2 = (estado == "av") ? (2.0 / 25.0) : (5.0 / 25.0)
    # Solo para verificacion: permite forzar el coeficiente de intercambio intra-capa
    # sin tocar la fisica por defecto. El valor por defecto (nothing) no cambia nada.
    coeff_k2_override === nothing || (coeff_k2 = coeff_k2_override)

    println("=== Silicon ROHF Optimización: $estado (Z=14) ===")
    
    # --- Physical & Grid Parameters ---
    N_elems = 100
    Z = 14.0
    
    ws = cached_init_scf_workspace(R_max, N_elems, Val(7), Z; γ=2.5)
    
    basis = ws.basis
    n = basis.num_splines
    active_s = 2:(n-1)  
    active_p = 3:(n-1)  
    
    # --- Orbital Configuration: 1s², 2s², 2p⁶, 3s², 3p² ---
    orbitals = [
        Orbital(1, 0, 2.0), # 1s (Closed)
        Orbital(2, 0, 2.0), # 2s (Closed)
        Orbital(2, 1, 6.0), # 2p (Closed Ne Core)
        Orbital(3, 0, 2.0), # 3s (Closed)
        Orbital(3, 1, 2.0)  # 3p (Open Valence, fractional w=2)
    ]
    
    H_core_s = ws.T + ws.V
    H_core_p = ws.T + ws.V + ws.R_inv2
    
    # --- Initial Guess Diagonalization ---
    evals_s, evecs_s = eigen(Symmetric(H_core_s[active_s, active_s]), ws.S[active_s, active_s])
    orbitals[1].coeffs = zeros(Float64, n); orbitals[1].coeffs[active_s] = evecs_s[:, 1]
    orbitals[2].coeffs = zeros(Float64, n); orbitals[2].coeffs[active_s] = evecs_s[:, 2]
    orbitals[4].coeffs = zeros(Float64, n); orbitals[4].coeffs[active_s] = evecs_s[:, 3]
    
    evals_p, evecs_p = eigen(Symmetric(H_core_p[active_p, active_p]), ws.S[active_p, active_p])
    orbitals[3].coeffs = zeros(Float64, n); orbitals[3].coeffs[active_p] = evecs_p[:, 1]
    orbitals[5].coeffs = zeros(Float64, n); orbitals[5].coeffs[active_p] = evecs_p[:, 2]
    
    for orb in orbitals
        orb.coeffs ./= sqrt(dot(orb.coeffs, ws.S * orb.coeffs))
    end

    # --- Pre-allocate Temporary SCF Matrices ---
    J_core         = zeros(Float64, n, n)
    K_core_s       = zeros(Float64, n, n)
    K_core_p       = zeros(Float64, n, n)
    
    J_3p_spherical = zeros(Float64, n, n)
    K_3p_on_s      = zeros(Float64, n, n)
    K_3p_on_p      = zeros(Float64, n, n)
    
    J_3p_intra     = zeros(Float64, n, n)
    K_3p_intra     = zeros(Float64, n, n)
    
    F_2p           = zeros(Float64, n, n)
    F_3p           = zeros(Float64, n, n)

    # --- SCF Loop Setup ---
    E_old = 0.0
    MIXING = 0.0
    diis = DIISState(; thresh = diis_thresh)
    diis_on = false
    # (iter, E_total, |dE|, residual proyectado, conmutador crudo) por iteracion:
    # es la curva que va a la figura de convergencia de resultados.tex.
    trace = NTuple{5,Float64}[]
    
    println("Starting ROHF SCF cycle...")
    if verbose
        @printf("%-4s | %-14s | %-10s | %-8s | %-10s | %s\n",
                "Iter", "E_total (Ha)", "Delta E", "Time (s)", "|FDS-SDF|", "modo")
        println("-"^78)
    end

    for iter in 1:max_iter
        t0 = time()

        # ==========================================
        # PHASE 1: Expanded Closed Core (1s + 2s + 2p + 3s)
        # ==========================================
        core_orbs = [orbitals[1], orbitals[2], orbitals[3], orbitals[4]]
        
        build_total_J_matrix!(ws, core_orbs)
        J_core .= ws.J
        
        assemble_K_matrix!(ws, ws.K_mats[0], 0, core_orbs)
        K_core_s .= ws.K_mats[0]
        
        assemble_K_matrix!(ws, ws.K_mats[1], 1, core_orbs)
        K_core_p .= ws.K_mats[1]

        # ==========================================
        # PHASE 2: Open 3p Shell Contributions
        # ==========================================
        
        # What the closed shells feel from the 3p shell
        build_specific_J_matrix!(ws, J_3p_spherical, orbitals[5], 2.0)
        
        assemble_K_matrix!(ws, ws.K_mats[0], 0, [orbitals[5]]) 
        K_3p_on_s .= ws.K_mats[0]
        
        assemble_K_matrix!(ws, ws.K_mats[1], 1, [orbitals[5]]) 
        K_3p_on_p .= ws.K_mats[1]

        # What the 3p electron feels from its own shell (Term-Dependent 3P)
        build_specific_J_matrix!(ws, J_3p_intra, orbitals[5], 1.0)
        
        # 3p Exchange Quadrupole (k=2, Coeff = 2/25 for av, -5/25 for 3P)
        build_specific_K_matrix!(ws, K_3p_intra, orbitals[5], 2, coeff_k2)

        # ==========================================
        # PHASE 3: Assemble Shell-Dependent Fock Matrices
        # ==========================================
        
        F_s = ws.F_mats[0]
        
        # s-block Operator (1s, 2s, 3s)
        F_s .= H_core_s .+ J_core .+ J_3p_spherical .- K_core_s .- K_3p_on_s
        
        # Core 2p Operator (Feels spherical 3p average)
        F_2p .= H_core_p .+ J_core .+ J_3p_spherical .- K_core_p .- K_3p_on_p
        
        # Valence 3p Operator (Feels intra-shell Fermi hole)
        F_3p .= H_core_p .+ J_core .+ J_3p_intra .- K_core_p .- K_3p_intra
        
        # ==========================================
        # PHASE 3b: C-DIIS por bloques (ver rohf_diis.jl)
        # ==========================================
        # Tres canales. F_2p y F_3p viven sobre el MISMO espacio activo l=1 pero son
        # operadores distintos, asi que cada uno lleva su propia densidad: la del 2p
        # cerrado (occ 6) y la del 3p abierto (occ 2).
        F_use = Matrix{Float64}[F_s, F_2p, F_3p]
        if use_diis
            D_s   = block_density(orbitals, (1, 2, 4), active_s)
            D_2p  = block_density(orbitals, (3,), active_p)
            D_3p  = block_density(orbitals, (5,), active_p)
            # Los dos canales l=1 comparten los mismos ocupados: es exactamente el
            # acoplamiento 2p-3p que hay que proyectar fuera del residual.
            C_s   = occ_block(orbitals, (1, 2, 4), active_s)
            C_p   = occ_block(orbitals, (3, 5), active_p)
            F_use, diis_on = diis_update!(diis, iter, Matrix{Float64}[F_s, F_2p, F_3p],
                                          Matrix{Float64}[D_s, D_2p, D_3p],
                                          Matrix{Float64}[C_s, C_p, C_p],
                                          UnitRange{Int}[active_s, active_p, active_p], ws.S)
        end
        F_s_use, F_2p_use, F_3p_use = F_use[1], F_use[2], F_use[3]

        # ==========================================
        # PHASE 4: Diagonalization & Projection
        # ==========================================
        
        evals_fs, evecs_fs   = eigen(Symmetric(F_s_use[active_s, active_s]), ws.S[active_s, active_s])
        evals_f2p, evecs_f2p = eigen(Symmetric(F_2p_use[active_p, active_p]), ws.S[active_p, active_p])
        evals_f3p, evecs_f3p = eigen(Symmetric(F_3p_use[active_p, active_p]), ws.S[active_p, active_p])
        
        # Extract new raw vectors
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
        
        # Update Orbitals (With Linear Mixing)
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
        # PHASE 5: Total Energy Calculation
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
        push!(trace, (Float64(iter), E_total, delta, diis.resid, diis.resid_raw))
        elapsed = time() - t0
        
        if verbose
            @printf("%-4d | %14.8f | %10.2e | %8.4f | %10.2e | %s\n",
                    iter, E_total, delta, elapsed, diis.resid, diis.frozen ? "libre" : (diis_on ? "DIIS" : "-"))
        end
        
        # No se acepta la convergencia mientras DIIS este extrapolando con el residual
        # ESTANCADO: ahi |dE| < tol puede dispararse sobre un ciclo limite en vez de
        # sobre un punto fijo. Con el residual todavia bajando (stall == 0) el minimo
        # es autentico y se acepta igual que sin DIIS.
        if delta < tol && (!diis.on || diis.frozen || diis.stall == 0)
            if verbose 
                println("-"^78)
                println("Converged in $iter iterations.")
            end
            @printf("Final HF Energy (%s): %.6f Ha\n", estado, E_total)

            if !save
                return (E_total = E_total, iters = iter, converged = true,
                        switch_iter = diis.switch_iter, freeze_iter = diis.freeze_iter,
                        restarts = diis.restarts,
                        resid = diis.resid, trace = trace, orbitals = orbitals)
            end

            dense_grid = exp.(range(log(1e-8), log(R_max), length=10000))
            V_eff = compute_effective_central_potential(ws, orbitals, dense_grid, Z)
            P_3p = evaluate_orbital(ws.basis, orbitals[5].coeffs, dense_grid)

            filename = "silicon_rohf_results_$(estado)_R$(R_max).jld2"
            jldsave(filename;
                    orbitals = orbitals,
                    E_total = E_total,
                    R_max = R_max,
                    R_grid = dense_grid,
                    V_nuclear = ws.V,
                    V_eff = V_eff,
                    P_3p = P_3p,
                    num_splines = n,
                    active_s = active_s,
                    active_p = active_p
                    )
            println("Saved Term-Dependent ROHF data to $filename")
            println("===== END =====")
            return (E_total = E_total, iters = iter, converged = true,
                    switch_iter = diis.switch_iter, freeze_iter = diis.freeze_iter,
                        restarts = diis.restarts,
                    resid = diis.resid, trace = trace, orbitals = orbitals)
        end
        
        E_old = E_total
    end
    @warn "Silicio: no convergio en $max_iter iteraciones (tol = $tol)."
    return (E_total = E_old, iters = max_iter, converged = false,
            switch_iter = diis.switch_iter, freeze_iter = diis.freeze_iter,
                        restarts = diis.restarts,
            resid = diis.resid, trace = trace, orbitals = orbitals)
end

if abspath(PROGRAM_FILE) == @__FILE__
    solve_silicon_rohf(30.0)
end
