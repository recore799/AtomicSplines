using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))

using AtomicSplines
using JLD2
using LinearAlgebra
using Printf
using WignerSymbols
using StaticArrays

include("../../src/ci.jl")

# Compute the angular coefficient for < l1^2 LS | 1/r12 | l2^2 LS >
function interaction_coefficient(l1::Int, l2::Int, L::Int, k::Int)
    w6j = wigner6j(Float64, l1, l1, L, l2, l2, k)
    w3j = wigner3j(Float64, l1, k, l2, 0, 0, 0)
    
    if abs(w3j) < 1e-10 || abs(w6j) < 1e-10
        return 0.0
    end
    
    rme2 = (2*l1 + 1) * (2*l2 + 1) * w3j^2
    # Fase correcta de <l'^2 LS|g|l^2 LS>: (-1)^(l+l'+L). En la diagonal (l1 == l2)
    # se reduce a (-1)^L, que era lo que habia antes.
    return (-1)^(l1 + l2 + L) * w6j * rme2
end

function get_h_core(ws::SolverWorkspace, a::Orbital, b::Orbital, core_orbs::Vector{Orbital})
    if a.l != b.l
        return 0.0
    end
    
    S_ab = dot(a.coeffs, ws.S * b.coeffs)
    
    # Core Hamiltonian: T + V_nuc + Centrifugal
    H_core_mat = ws.T .+ ws.V .+ ((a.l * (a.l + 1)) / 2.0) .* ws.R_inv2
    
    val = dot(a.coeffs, H_core_mat * b.coeffs)
    
    # Add repulsion from each core orbital
    for c in core_orbs
        # Direct term J
        # R^0(a, c, b, c) * 2(2l_c+1)
        val += 2.0 * (2*c.l + 1) * get_cached_Rk!(ws, a, c, b, c, 0)
        
        # Exchange term K
        # - sum_k c_k R^k(a, c, c, b)
        for k in abs(a.l - c.l):(a.l + c.l)
            w3j = wigner3j(Float64, a.l, k, c.l, 0, 0, 0)
            if abs(w3j) > 1e-10
                coeff = (2*c.l + 1) * w3j^2
                val -= coeff * get_cached_Rk!(ws, a, c, c, b, k)
            end
        end
    end
    
    return val
end

const NP2_VALENCE_N = Dict("Carbon" => 2, "Silicon" => 3, "Germanium" => 4, "Tin" => 5)

"""
    compute_zeta(dense_grid, V_eff, P_nl)

Parametro radial de espin-orbita, zeta = (alpha^2/2) <P_nl| (1/r) dV/dr |P_nl>.
Misma rutina que usan `*_fine_structure.jl` y los scripts de calibracion de V_pol, para
que el CI no dependa de un literal escrito a mano.
"""
function compute_zeta(dense_grid::Vector{Float64}, V_eff::Vector{Float64}, P_nl::Vector{Float64})
    N = length(dense_grid)
    dV_dr = zeros(Float64, N)
    integrand = zeros(Float64, N)

    alpha = 1.0 / 137.035999

    dV_dr[1] = (V_eff[2] - V_eff[1]) / (dense_grid[2] - dense_grid[1])
    dV_dr[N] = (V_eff[N] - V_eff[N-1]) / (dense_grid[N] - dense_grid[N-1])
    for i in 2:(N-1)
        dV_dr[i] = (V_eff[i+1] - V_eff[i-1]) / (dense_grid[i+1] - dense_grid[i-1])
    end

    for i in 1:N
        r = dense_grid[i]
        integrand[i] = r > 1e-12 ? (P_nl[i]^2) * (1.0 / r) * dV_dr[i] : 0.0
    end

    zeta_val = 0.0
    for i in 1:(N-1)
        dr = dense_grid[i+1] - dense_grid[i]
        zeta_val += 0.5 * (integrand[i] + integrand[i+1]) * dr
    end

    return (alpha^2 / 2.0) * zeta_val
end

"""
    s_orthonormalize!(orbs, ws)

Gram-Schmidt modificado en la metrica de traslape S. Hace falta porque el orbital SCF de
valencia se inserta en un juego de autovectores del Fock de core congelado: esos virtuales
son ortogonales entre si y al "np" congelado, pero no al np del SCF completo. El CI
diagonaliza asumiendo <i|S|j> = delta_ij, asi que la base tiene que serlo de verdad.
"""
function s_orthonormalize!(orbs::Vector{Orbital}, ws::SolverWorkspace)
    for i in eachindex(orbs)
        ci = orbs[i].coeffs
        for j in 1:(i-1)
            cj = orbs[j].coeffs
            ci .-= dot(cj, ws.S * ci) .* cj
        end
        nrm = sqrt(dot(ci, ws.S * ci))
        if nrm < 1e-8
            error("Espacio activo linealmente dependiente en l=$(orbs[i].l), i=$i (norma $nrm).")
        end
        ci ./= nrm
    end
    return orbs
end

"""
    pick_reference_root(evals, evecs, idx_ref)

(energia, |coeficiente|, indice) de la raiz con mayor peso sobre el CSF de referencia.
No se toma la mas baja del bloque: en cuanto el espacio activo crece, o V_pol mueve las
energias, la raiz mas baja de ^1S puede ser un s^2 o un Rydberg y el bloque J=0 apuntaria
al estado equivocado sin avisar. La fase se fija con abs() porque eigen() no la garantiza.
"""
function pick_reference_root(evals::Vector{Float64}, evecs::Matrix{Float64}, idx_ref::Int)
    j = argmax(abs.(@view evecs[idx_ref, :]))
    return evals[j], abs(evecs[idx_ref, j]), j
end


function run_toy_ci(element_name::String, result_file::String, zeta_np::Float64)
    println("==================================================")
    println(" Toy CI for $(element_name) (Valence correlation only) ")
    println("==================================================")
    
    data = load(joinpath(@__DIR__, result_file))
    R_max = data["R_max"]
    orbitals = data["orbitals"]
    
    alpha_d = get(data, "alpha_d", 0.0)
    r_c = get(data, "r_c", 1.0)
    
    local ws
    local Z::Float64
    if element_name == "Carbon"
        Z = 6.0
        ws = cached_init_scf_workspace(R_max, 100, Val(7), Z; γ=2.5, alpha_d=alpha_d, r_c=r_c)
    elseif element_name == "Silicon"
        Z = 14.0
        ws = cached_init_scf_workspace(R_max, 100, Val(7), Z; γ=2.5, alpha_d=alpha_d, r_c=r_c)
    elseif element_name == "Germanium"
        Z = 32.0
        ws = cached_init_scf_workspace(R_max, 300, Val(8), Z; γ=3.0, alpha_d=alpha_d, r_c=r_c)
    elseif element_name == "Tin"
        Z = 50.0
        ws = cached_init_scf_workspace(R_max, 500, Val(8), Z; γ=4.0, alpha_d=alpha_d, r_c=r_c)
    else
        error("Element not supported: $element_name")
    end
    
    val_orb = orbitals[end]
    core_orbs = orbitals[1:end-1]
    
    np_orb = val_orb
    
    println("Building Frozen Core Fock Operator to extract virtuals...")
    build_total_J_matrix!(ws, core_orbs)
    J_core = copy(ws.J)
    
    H_core_s = ws.T .+ ws.V .+ ((0.0) / 2.0) .* ws.R_inv2 .+ J_core
    H_core_p = ws.T .+ ws.V .+ ((2.0) / 2.0) .* ws.R_inv2 .+ J_core
    H_core_d = ws.T .+ ws.V .+ ((6.0) / 2.0) .* ws.R_inv2 .+ J_core
    
    assemble_K_matrix!(ws, ws.K_mats[0], 0, core_orbs)
    F_s = H_core_s .- ws.K_mats[0]
    evals_s, evecs_s = eigen(Symmetric(F_s), Symmetric(ws.S))
    
    assemble_K_matrix!(ws, ws.K_mats[1], 1, core_orbs)
    F_p = H_core_p .- ws.K_mats[1]
    evals_p, evecs_p = eigen(Symmetric(F_p), Symmetric(ws.S))
    
    assemble_K_matrix!(ws, ws.K_mats[2], 2, core_orbs)
    F_d = H_core_d .- ws.K_mats[2]
    evals_d, evecs_d = eigen(Symmetric(F_d), Symmetric(ws.S))
    
    # Active space: the actual np orbital + virtual p + virtual d + virtual s
    active_s = extract_virtuals(evals_s, evecs_s, 0, count(o->o.l==0, core_orbs), 2, 1:ws.basis.num_splines, ws.basis.num_splines, ws, Z)
    active_p = extract_virtuals(evals_p, evecs_p, 1, count(o->o.l==1, core_orbs), 3, 1:ws.basis.num_splines, ws.basis.num_splines, ws, Z; offset=0)
    active_d = extract_virtuals(evals_d, evecs_d, 2, count(o->o.l==2, core_orbs), 2, 1:ws.basis.num_splines, ws.basis.num_splines, ws, Z)
    
    # El orbital SCF de valencia no es autovector del Fock de core congelado, asi que al
    # insertarlo la base deja de ser ortonormal. Se copia (para no mutar el orbital cargado
    # del .jld2) y se reortonormaliza en la metrica S.
    active_p[1] = deepcopy(np_orb)
    s_orthonormalize!(active_s, ws)
    s_orthonormalize!(active_p, ws)
    s_orthonormalize!(active_d, ws)
    
    csfs = vcat(active_s, active_p, active_d)
    
    println("Active Space CSFs (n l^2):")
    for (i, orb) in enumerate(csfs)
        println("  $i: $(orb.n)$(orb.l == 0 ? "s" : (orb.l == 1 ? "p" : "d"))^2")
    end
    
    # Build CI Matrix for 3P (L=1, S=1)
    csfs_3P = filter(o -> o.l > 0, csfs)
    N_3P = length(csfs_3P)
    H_3P = zeros(Float64, N_3P, N_3P)
    
    for i in 1:N_3P
        orb_a = csfs_3P[i]
        H_3P[i, i] = 2.0 * get_h_core(ws, orb_a, orb_a, core_orbs)
        
        for k in 0:2:(2*orb_a.l)
            c_k = interaction_coefficient(orb_a.l, orb_a.l, 1, k)
            H_3P[i, i] += c_k * get_cached_Rk!(ws, orb_a, orb_a, orb_a, orb_a, k)
        end
        
        for j in (i+1):N_3P
            orb_b = csfs_3P[j]
            val = 0.0
            for k in abs(orb_a.l - orb_b.l):2:(orb_a.l + orb_b.l)
                c_k = interaction_coefficient(orb_a.l, orb_b.l, 1, k)
                val += c_k * get_cached_Rk!(ws, orb_a, orb_a, orb_b, orb_b, k)
            end
            H_3P[i, j] = val
            H_3P[j, i] = val
        end
    end
    
    # Build CI Matrix for 1D (L=2, S=0)
    csfs_1D = filter(o -> o.l > 0, csfs)
    N_1D = length(csfs_1D)
    H_1D = zeros(Float64, N_1D, N_1D)
    
    for i in 1:N_1D
        orb_a = csfs_1D[i]
        H_1D[i, i] = 2.0 * get_h_core(ws, orb_a, orb_a, core_orbs)
        
        for k in 0:2:(2*orb_a.l)
            c_k = interaction_coefficient(orb_a.l, orb_a.l, 2, k)
            H_1D[i, i] += c_k * get_cached_Rk!(ws, orb_a, orb_a, orb_a, orb_a, k)
        end
        
        for j in (i+1):N_1D
            orb_b = csfs_1D[j]
            val = 0.0
            for k in abs(orb_a.l - orb_b.l):2:(orb_a.l + orb_b.l)
                c_k = interaction_coefficient(orb_a.l, orb_b.l, 2, k)
                val += c_k * get_cached_Rk!(ws, orb_a, orb_a, orb_b, orb_b, k)
            end
            H_1D[i, j] = val
            H_1D[j, i] = val
        end
    end
    
    # Build CI Matrix for 1S (L=0, S=0)
    csfs_1S = csfs
    N_1S = length(csfs_1S)
    H_1S = zeros(Float64, N_1S, N_1S)
    
    for i in 1:N_1S
        orb_a = csfs_1S[i]
        H_1S[i, i] = 2.0 * get_h_core(ws, orb_a, orb_a, core_orbs)
        
        for k in 0:2:(2*orb_a.l)
            c_k = interaction_coefficient(orb_a.l, orb_a.l, 0, k)
            H_1S[i, i] += c_k * get_cached_Rk!(ws, orb_a, orb_a, orb_a, orb_a, k)
        end
        
        for j in (i+1):N_1S
            orb_b = csfs_1S[j]
            val = 0.0
            for k in abs(orb_a.l - orb_b.l):2:(orb_a.l + orb_b.l)
                c_k = interaction_coefficient(orb_a.l, orb_b.l, 0, k)
                val += c_k * get_cached_Rk!(ws, orb_a, orb_a, orb_b, orb_b, k)
            end
            H_1S[i, j] = val
            H_1S[j, i] = val
        end
    end
    
    # --- Diagonalizacion de los bloques LS ---
    println("H_3P diagonal:")
    println(diag(H_3P))
    println("H_1S diagonal:")
    println(diag(H_1S))

    e_3P, v_3P = eigen(Symmetric(H_3P))
    e_1D, v_1D = eigen(Symmetric(H_1D))
    e_1S, v_1S = eigen(Symmetric(H_1S))

    # Indice del CSF de referencia np^2 dentro de cada bloque
    idx_np_3P = findfirst(o -> o.l == 1 && o.n == np_orb.n, csfs_3P)
    idx_np_1D = findfirst(o -> o.l == 1 && o.n == np_orb.n, csfs_1D)
    idx_np_1S = findfirst(o -> o.l == 1 && o.n == np_orb.n, csfs_1S)

    E_av_CI, C_3P, root_3P = pick_reference_root(e_3P, v_3P, idx_np_3P)
    E_1D_CI, C_1D, root_1D = pick_reference_root(e_1D, v_1D, idx_np_1D)
    E_1S_CI, C_1S, root_1S = pick_reference_root(e_1S, v_1S, idx_np_1S)

    for (nm, r) in (("3P", root_3P), ("1D", root_1D), ("1S", root_1S))
        if r != 1
            @warn "La raiz de $(nm) con mayor peso np^2 es la #$(r), no la mas baja del bloque."
        end
    end

    # --- Test de consistencia angular ---
    # Sobre el CSF de referencia los coeficientes f_k son exactos y no dependen del modelo:
    #   E(^1D) - E(^3P) = (1/25 + 5/25) F^2 = 0.24 F^2
    #   E(^1S) - E(^3P) = (10/25 + 5/25) F^2 = 0.60 F^2
    # Si esto falla, el problema esta en compute_Rk o en interaction_coefficient, y todo lo
    # que sigue (incluido el bloque J=0 de Breit-Pauli) hereda el error.
    ref_orb = csfs_3P[idx_np_3P]
    F2_ref = get_cached_Rk!(ws, ref_orb, ref_orb, ref_orb, ref_orb, 2)
    d_3P = H_3P[idx_np_3P, idx_np_3P]
    d_1D = H_1D[idx_np_1D, idx_np_1D]
    d_1S = H_1S[idx_np_1S, idx_np_1S]
    r_1D = (d_1D - d_3P) / F2_ref
    r_1S = (d_1S - d_3P) / F2_ref
    angular_ok = abs(r_1D - 0.24) <= 1e-8 && abs(r_1S - 0.60) <= 1e-8

    println("\n--- Test de consistencia angular (CSF de referencia) ---")
    @printf("  F^2(%dp,%dp)         = %.8f Ha (%.2f cm^-1)\n",
            np_orb.n, np_orb.n, F2_ref, F2_ref * 219474.63)
    @printf("  [E(1D)-E(3P)]/F^2   = %.10f   exacto 0.24   delta %.2e\n", r_1D, r_1D - 0.24)
    @printf("  [E(1S)-E(3P)]/F^2   = %.10f   exacto 0.60   delta %.2e\n", r_1S, r_1S - 0.60)
    if !angular_ok
        @warn "Los coeficientes angulares NO se reproducen. Revisar compute_Rk / interaction_coefficient ANTES de usar estos numeros."
    end

    println("\n--- CI Results ---")
    @printf("  zeta_np usado       = %.8f Ha (%.3f cm^-1)\n", zeta_np, zeta_np * 219474.63)
    println("E(3P) = ", E_av_CI)
    println("E_corr(3P) = ", E_av_CI - d_3P)
    @printf("  peso np^2: C(3P) = %.6f  C(1D) = %.6f  C(1S) = %.6f\n", C_3P, C_1D, C_1S)
    println("E(1D) = ", E_1D_CI)
    println("E(1S) = ", E_1S_CI)
    println("Splitting 1D - 3P = ", (E_1D_CI - E_av_CI) * 219474.63, " cm-1")
    println("Splitting 1S - 3P = ", (E_1S_CI - E_av_CI) * 219474.63, " cm-1")
    
    # Breit-Pauli matrices
    # J = 0 Manifold (3P0, 1S0)
    H_BP0 = zeros(Float64, 2, 2)
    H_BP0[1, 1] = E_av_CI - zeta_np * C_3P^2
    H_BP0[2, 2] = E_1S_CI
    H_BP0[1, 2] = sqrt(2.0) * zeta_np * C_3P * C_1S
    H_BP0[2, 1] = H_BP0[1, 2]
    evals_BP0, evecs_BP0 = eigen(Symmetric(H_BP0))
    E_av_0_CI = evals_BP0[1]
    E_1S_0_CI = evals_BP0[2]
    
    # J = 1 Manifold (3P1)
    E_av_1_CI = E_av_CI - 0.5 * zeta_np * C_3P^2
    
    # J = 2 Manifold (3P2, 1D2)
    H_BP2 = zeros(Float64, 2, 2)
    H_BP2[1, 1] = E_av_CI + 0.5 * zeta_np * C_3P^2
    H_BP2[2, 2] = E_1D_CI
    H_BP2[1, 2] = -(1.0/sqrt(2.0)) * zeta_np * C_3P * C_1D
    H_BP2[2, 1] = H_BP2[1, 2]
    evals_BP2, evecs_BP2 = eigen(Symmetric(H_BP2))
    E_av_2_CI = evals_BP2[1]
    c1_3P = evecs_BP2[1, 1]
    c1_1D = evecs_BP2[2, 1]
    g_eff = c1_3P^2 * 1.50 + c1_1D^2 * 1.00
    println("g_eff_CI (3P2) = ", g_eff)
    H_BP2_HF = zeros(Float64, 2, 2)
    H_BP2_HF[1, 1] = H_3P[idx_np_3P, idx_np_3P] + 0.5 * zeta_np
    H_BP2_HF[2, 2] = H_1D[idx_np_1D, idx_np_1D]
    H_BP2_HF[1, 2] = -(1.0/sqrt(2.0)) * zeta_np
    H_BP2_HF[2, 1] = H_BP2_HF[1, 2]
    evals_BP2_HF, evecs_BP2_HF = eigen(Symmetric(H_BP2_HF))
    g_eff_HF = evecs_BP2_HF[1, 1]^2 * 1.50 + evecs_BP2_HF[2, 1]^2 * 1.00
    println("g_eff_HF (3P2) = ", g_eff_HF)
    E_1D_2_CI = evals_BP2[2]
    
    au2cm = 219474.63
    rel_3P_0 = 0.0
    rel_3P_1 = (E_av_1_CI - E_av_0_CI) * au2cm
    rel_3P_2 = (E_av_2_CI - E_av_0_CI) * au2cm
    rel_1D_2 = (E_1D_2_CI - E_av_0_CI) * au2cm
    rel_1S_0 = (E_1S_0_CI - E_av_0_CI) * au2cm
    
    println("--- Corrected Multiplet Spectrum (cm-1) ---")
    println("3P_0 : ", rel_3P_0)
    println("3P_1 : ", rel_3P_1)
    println("3P_2 : ", rel_3P_2)
    println("1D_2 : ", rel_1D_2)
    println("1S_0 : ", rel_1S_0)
    println("-------------------------------------------")

    return (levels     = [rel_3P_0, rel_3P_1, rel_3P_2, rel_1D_2, rel_1S_0],
            E_3P       = E_av_CI,
            E_1D       = E_1D_CI,
            E_1S       = E_1S_CI,
            E_corr     = E_av_CI - d_3P,
            F2         = F2_ref,
            zeta       = zeta_np,
            C          = (C_3P, C_1D, C_1S),
            angular_ok = angular_ok,
            g_eff      = g_eff,
            g_eff_HF   = g_eff_HF)
end

"""
    run_toy_ci(element_name, result_file)

Igual que el metodo de tres argumentos, pero calculando zeta_np del propio archivo SCF
(R_grid, V_eff y P_np) en vez de recibirlo como literal. Esta es la forma correcta de
llamarlo al barrer alpha_d: V_pol relaja el orbital, baja <1/r^3> y con ello zeta. Con un
zeta literal el barrido no mueve la estructura fina y parece que V_pol no hace nada.
"""
function run_toy_ci(element_name::String, result_file::String)
    haskey(NP2_VALENCE_N, element_name) || error("Element not supported: $element_name")
    n_val = NP2_VALENCE_N[element_name]
    data = load(joinpath(@__DIR__, result_file))
    key = "P_$(n_val)p"
    haskey(data, key) || error("$(result_file) no contiene la clave $(key); regenera el .jld2 con el script ROHF correspondiente.")
    zeta_np = compute_zeta(data["R_grid"], data["V_eff"], data[key])
    return run_toy_ci(element_name, result_file, zeta_np)
end


# Execute tests
# NOTE: The CI is best run on Configuration Average ('av') orbitals so that the basis is
# not biased towards the 3P ground state.
# zeta se calcula del .jld2; ya no se pasa como literal.
if abspath(PROGRAM_FILE) == @__FILE__
    run_toy_ci("Carbon",    "carbon_rohf_results_3P_R30.0.jld2")
    run_toy_ci("Silicon",   "silicon_rohf_results_3P_R30.0.jld2")
    run_toy_ci("Germanium", "germanium_rohf_results_3P_R30.0.jld2")
    run_toy_ci("Tin",       "tin_rohf_results_3P_R30.0.jld2")
end