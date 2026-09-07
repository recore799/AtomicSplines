# =============================================================================
#  CI de pareja de valencia para configuraciones np^2, general en (n l)
#
#  Sustituye al CI de "pares equivalentes" de np2_toy_ci.jl, que solo admitia
#  CSFs del tipo (nl)^2 y por eso recuperaba ~1-5% de la correlacion de valencia.
#  Aqui la base son todos los CSFs |(n1 l1, n2 l2) L S> antisimetrizados que
#  respetan paridad par, la regla del triangulo y Pauli.
#
#  ADVERTENCIA DE COSTO: compute_Rk resuelve una ecuacion de Poisson por integral,
#  sin lista transformada de integrales bi-electronicas. El numero de R^k distintos
#  crece como (N_orb)^4 / simetria, asi que el "pseudo-espectro completo" de splines
#  NO es alcanzable. El uso previsto es un ESTUDIO DE CONVERGENCIA: subir n_per_l
#  hasta que E_corr se estabilice, y reportar esa curva.
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))

using AtomicSplines
using JLD2
using LinearAlgebra
using Printf
using WignerSymbols

# get_h_core, compute_zeta, s_orthonormalize!, pick_reference_root,
# NP2_VALENCE_N, interaction_coefficient. El bloque de ejecucion de ese archivo
# esta protegido por abspath(PROGRAM_FILE) == @__FILE__, asi que no corre solo.
include("np2_toy_ci.jl")

const AU2CM = 219474.63
const LCHAR = ['s', 'p', 'd', 'f', 'g', 'h', 'i']

orb_label(o::Orbital) = string(o.n) * string(LCHAR[o.l + 1])

# -----------------------------------------------------------------------------
#  CSF: pareja antisimetrizada |(i, j) L S> con i <= j sobre el pool de orbitales
# -----------------------------------------------------------------------------
struct PairCSF
    i::Int
    j::Int
    L::Int
    S::Int
end

csf_label(c::PairCSF, pool) = orb_label(pool[c.i]) * orb_label(pool[c.j])

is_equivalent(c::PairCSF) = c.i == c.j

# Norma del par antisimetrizado: 1/2 para electrones equivalentes, 1/sqrt(2) si no.
norm_fac(c::PairCSF) = is_equivalent(c) ? 0.5 : (1.0 / sqrt(2.0))

# P_12 |ab;LS> = -(-1)^(la+lb+L+S) |ba;LS>, de donde la fase de intercambio.
xchg_phase(pool, c::PairCSF) = (-1)^(pool[c.i].l + pool[c.j].l + c.L + c.S)

"""
    generate_pair_csfs(pool, L, S; parity)

Todos los CSFs de dos electrones con momento angular total (L,S) y la paridad dada
(`parity = mod(l1+l2, 2)`; 0 para np^2). Para electrones equivalentes (i == j) solo
sobrevive L+S par, que es el principio de Pauli en esta base acoplada.
"""
function generate_pair_csfs(pool::Vector{Orbital}, L::Int, S::Int; parity::Int = 0)
    out = PairCSF[]
    n = length(pool)
    for i in 1:n, j in i:n
        la, lb = pool[i].l, pool[j].l
        mod(la + lb, 2) == parity || continue
        (abs(la - lb) <= L <= la + lb) || continue
        if i == j && !iseven(L + S)
            continue
        end
        push!(out, PairCSF(i, j, L, S))
    end
    return out
end

# -----------------------------------------------------------------------------
#  Elementos de matriz
# -----------------------------------------------------------------------------

"""
    h_direct(ws, pool, core, a, b, c, d)

<ab| h(1) + h(2) |cd> sin antisimetrizar, con h el hamiltoniano de core congelado
(cinetica + nuclear + centrifugo + J_core - K_core). Los orbitales del pool son
S-ortonormales dentro de cada l, de modo que <b|d> = delta_bd.
"""
function h_direct(ws, pool, core, a::Int, b::Int, c::Int, d::Int)
    (pool[a].l == pool[c].l && pool[b].l == pool[d].l) || return 0.0
    v = 0.0
    if b == d
        v += get_h_core(ws, pool[a], pool[c], core)
    end
    if a == c
        v += get_h_core(ws, pool[b], pool[d], core)
    end
    return v
end

"""
    g_direct(ws, pool, a, b, c, d, L)

Termino directo <ab| 1/r12 |cd> acoplado a L:

    sum_k (-1)^(lb+lc+L) { la lb L ; ld lc k } <la||C^k||lc> <lb||C^k||ld> R^k(ab,cd)

con <l||C^k||l'> = (-1)^l sqrt((2l+1)(2l'+1)) (l k l'; 0 0 0).
El intercambio se obtiene llamando con (c,d) permutados, ver pair_matrix_element.
"""
function g_direct(ws, pool, a::Int, b::Int, c::Int, d::Int, L::Int)
    la, lb = pool[a].l, pool[b].l
    lc, ld = pool[c].l, pool[d].l
    kmin = max(abs(la - lc), abs(lb - ld))
    kmax = min(la + lc, lb + ld)
    tot = 0.0
    for k in kmin:kmax
        w3a = wigner3j(Float64, la, k, lc, 0, 0, 0)
        abs(w3a) < 1e-12 && continue
        w3b = wigner3j(Float64, lb, k, ld, 0, 0, 0)
        abs(w3b) < 1e-12 && continue
        w6 = wigner6j(Float64, la, lb, L, ld, lc, k)
        abs(w6) < 1e-12 && continue
        rme = (-1)^(la + lb) *
              sqrt((2la + 1) * (2lc + 1) * (2lb + 1) * (2ld + 1)) * w3a * w3b
        tot += (-1)^(lb + lc + L) * w6 * rme *
               get_cached_Rk!(ws, pool[a], pool[b], pool[c], pool[d], k)
    end
    return tot
end

"""
    pair_matrix_element(ws, pool, core, c1, c2)

<c1| H |c2> entre CSFs antisimetrizados y normalizados:

    2 N_ab N_cd [ <ab|H|cd> + (-1)^(lc+ld+L+S) <ab|H|dc> ]

La fase del bra y la del ket coinciden porque toda la base tiene la misma paridad.
"""
function pair_matrix_element(ws, pool, core, c1::PairCSF, c2::PairCSF)
    (c1.L == c2.L && c1.S == c2.S) || return 0.0
    a, b = c1.i, c1.j
    c, d = c2.i, c2.j
    L = c1.L
    ph = xchg_phase(pool, c2)
    pref = 2.0 * norm_fac(c1) * norm_fac(c2)
    one_body = h_direct(ws, pool, core, a, b, c, d) +
               ph * h_direct(ws, pool, core, a, b, d, c)
    two_body = g_direct(ws, pool, a, b, c, d, L) +
               ph * g_direct(ws, pool, a, b, d, c, L)
    return pref * (one_body + two_body)
end

function build_pair_hamiltonian(ws, pool, core, csfs::Vector{PairCSF}; verbose::Bool = false)
    n = length(csfs)
    H = zeros(Float64, n, n)
    for i in 1:n
        for j in i:n
            v = pair_matrix_element(ws, pool, core, csfs[i], csfs[j])
            H[i, j] = v
            H[j, i] = v
        end
        if verbose && (i % 25 == 0 || i == n)
            @printf("    fila %d/%d  (R^k en cache: %d)\n", i, n, length(ws.rk_cache))
        end
    end
    return H
end

# -----------------------------------------------------------------------------
#  Pool de orbitales
# -----------------------------------------------------------------------------
"""
    build_orbital_pool(ws, core_orbs, val_orb, Z; n_per_l)

Diagonaliza el Fock de core congelado (V^{N-2}) por cada l y toma los primeros
`n_per_l[l]` virtuales fisicos. En el canal de valencia el primero se sustituye por
el orbital SCF real, y despues cada bloque se reortonormaliza en la metrica S
(los virtuales congelados son ortogonales al "np" congelado, no a este).
"""
function build_orbital_pool(ws, core_orbs::Vector{Orbital}, val_orb::Orbital, Z::Float64;
                            n_per_l::Dict{Int,Int} = Dict(0 => 6, 1 => 6, 2 => 6, 3 => 4))
    build_total_J_matrix!(ws, core_orbs)
    J_core = copy(ws.J)

    pool = Orbital[]
    for l in sort(collect(keys(n_per_l)))
        nvirt = n_per_l[l]
        nvirt > 0 || continue
        l <= 3 || error("assemble_K_matrix! solo cubre l <= 3; pediste l = $l.")

        H_l = ws.T .+ ws.V .+ ((l * (l + 1)) / 2.0) .* ws.R_inv2 .+ J_core
        assemble_K_matrix!(ws, ws.K_mats[l], l, core_orbs)
        F_l = H_l .- ws.K_mats[l]
        ev, vec = eigen(Symmetric(F_l), Symmetric(ws.S))

        block = extract_virtuals(ev, vec, l, count(o -> o.l == l, core_orbs), nvirt,
                                 1:ws.basis.num_splines, ws.basis.num_splines, ws, Z)
        length(block) == nvirt ||
            @warn "l=$l: se pidieron $nvirt virtuales y solo hay $(length(block)) fisicos."

        if l == val_orb.l && !isempty(block)
            block[1] = deepcopy(val_orb)
        end
        s_orthonormalize!(block, ws)
        append!(pool, block)
    end
    return pool
end

# -----------------------------------------------------------------------------
#  Autoverificacion
# -----------------------------------------------------------------------------
"""
    engine_selftest(ws, pool, core, val_idx)

Tres pruebas que no dependen de datos externos:

  T1  Hermiticidad de cada bloque (L,S), reevaluando el triangulo superior con los
      CSFs invertidos (la matriz de build_pair_hamiltonian ya sale simetrizada).
  T2  Sobre el subespacio de pares equivalentes (nl)^2 el motor general debe
      reproducir DIGITO A DIGITO al motor legado de np2_toy_ci.jl, que usa
      interaction_coefficient. Valida fases, normas y 6j de ESE subespacio.
  T3  Coeficientes exactos de p^2 sobre el CSF de referencia:
      [E(1D)-E(3P)]/F^2 = 0.24 y [E(1S)-E(3P)]/F^2 = 0.60.
  T4  Pares NO equivalentes (ns, n'l') contra Condon-Shortley. Es la unica prueba
      que toca lo que T2 no puede tocar, que es justo lo que el motor agrega.

Devuelve true si todo pasa.
"""
function engine_selftest(ws, pool, core, val_idx::Int; tol::Float64 = 1e-10)
    ok = true
    println("\n===== Autoverificacion del motor CI =====")

    # --- T1 y T2 ---
    for (L, S, name) in ((1, 1, "3P"), (2, 0, "1D"), (0, 0, "1S"))
        csfs = generate_pair_csfs(pool, L, S)
        isempty(csfs) && continue
        H = build_pair_hamiltonian(ws, pool, core, csfs)

        # build_pair_hamiltonian llena H[p,q] y H[q,p] con la MISMA llamada, de modo
        # que maximum(abs.(H .- H')) sobre esa matriz vale cero por construccion y no
        # prueba nada. La hermiticidad hay que medirla reevaluando el triangulo
        # superior con los dos CSFs en orden invertido.
        asym = 0.0
        for p in eachindex(csfs), q in (p + 1):length(csfs)
            asym = max(asym,
                       abs(H[p, q] - pair_matrix_element(ws, pool, core, csfs[q], csfs[p])))
        end
        if asym > tol
            @warn "T1 $name: hermiticidad rota, max|H-H'| = $asym"
            ok = false
        end

        # T2: comparar solo entre CSFs equivalentes
        eq = [k for k in eachindex(csfs) if is_equivalent(csfs[k])]
        worst = 0.0
        for p in eq, q in eq
            ca, cb = csfs[p], csfs[q]
            oa, ob = pool[ca.i], pool[cb.i]
            ref = 0.0
            if ca.i == cb.i
                ref += 2.0 * get_h_core(ws, oa, oa, core)
                for k in 0:2:(2 * oa.l)
                    ref += interaction_coefficient(oa.l, oa.l, L, k) *
                           get_cached_Rk!(ws, oa, oa, oa, oa, k)
                end
            else
                for k in abs(oa.l - ob.l):2:(oa.l + ob.l)
                    ref += interaction_coefficient(oa.l, ob.l, L, k) *
                           get_cached_Rk!(ws, oa, oa, ob, ob, k)
                end
            end
            worst = max(worst, abs(H[p, q] - ref))
        end
        status = worst <= tol ? "OK" : "FALLA"
        @printf("  T1/T2 %-3s : %3d CSFs (%2d equivalentes)  asim %.2e  vs legado %.2e  [%s]\n",
                name, length(csfs), length(eq), asym, worst, status)
        worst > tol && (ok = false)
    end

    # --- T3 ---
    diags = Dict{String,Float64}()
    for (L, S, name) in ((1, 1, "3P"), (2, 0, "1D"), (0, 0, "1S"))
        csfs = generate_pair_csfs(pool, L, S)
        idx = findfirst(c -> c.i == val_idx && c.j == val_idx, csfs)
        idx === nothing && error("El CSF de referencia no aparece en el bloque $name.")
        c = csfs[idx]
        diags[name] = pair_matrix_element(ws, pool, core, c, c)
    end
    vo = pool[val_idx]
    F2 = get_cached_Rk!(ws, vo, vo, vo, vo, 2)
    r1D = (diags["1D"] - diags["3P"]) / F2
    r1S = (diags["1S"] - diags["3P"]) / F2
    @printf("  T3 angular : F^2 = %.8f Ha | [1D-3P]/F^2 = %.10f (0.24, d=%.1e) | [1S-3P]/F^2 = %.10f (0.60, d=%.1e)\n",
            F2, r1D, r1D - 0.24, r1S, r1S - 0.60)
    if abs(r1D - 0.24) > 1e-8 || abs(r1S - 0.60) > 1e-8
        @warn "T3 FALLA: los coeficientes angulares exactos de p^2 no se reproducen."
        ok = false
    end

    # --- T4 ---
    # Pares NO equivalentes contra las formulas cerradas de Condon-Shortley. T2 solo
    # compara el subespacio (nl)^2, que es exactamente lo que ya hacia el motor legado:
    # sin T4 la parte nueva del motor -- normas 1/sqrt(2), fase de intercambio y 6j con
    # l_a != l_b -- se queda sin verificar. Para la configuracion (ns, n'l') sobreviven
    # solo el directo F^0 y el intercambio G^l', y los dos terminos valen
    #     E(^1L) = I_a + I_b + F^0(ab) + G^l'(ab)/(2l'+1)
    #     E(^3L) = I_a + I_b + F^0(ab) - G^l'(ab)/(2l'+1),   con L = l'.
    ia = findfirst(o -> o.l == 0, pool)
    if ia === nothing
        println("  T4 no-equiv: el pool no tiene orbitales s, prueba omitida.")
    else
        oa = pool[ia]
        Ia = get_h_core(ws, oa, oa, core)
        for lb in 0:maximum(o.l for o in pool)
            ib = findfirst(k -> k != ia && pool[k].l == lb, eachindex(pool))
            ib === nothing && continue
            ob = pool[ib]
            p, q = minmax(ia, ib)
            F0 = get_cached_Rk!(ws, oa, ob, oa, ob, 0)
            Gl = get_cached_Rk!(ws, oa, ob, ob, oa, lb)
            Ib = get_h_core(ws, ob, ob, core)
            worst = 0.0
            for (S, sgn) in ((0, 1.0), (1, -1.0))
                cs = generate_pair_csfs(pool, lb, S; parity = mod(lb, 2))
                r = findfirst(c -> c.i == p && c.j == q, cs)
                if r === nothing
                    @warn "T4: falta el CSF $(orb_label(oa))$(orb_label(ob)) en L=$lb, S=$S."
                    ok = false
                    continue
                end
                got = pair_matrix_element(ws, pool, core, cs[r], cs[r])
                worst = max(worst, abs(got - (Ia + Ib + F0 + sgn * Gl / (2 * lb + 1))))
            end
            status = worst <= tol ? "OK" : "FALLA"
            @printf("  T4 %-6s : F^0 = %.8f  G^%d = %.8f  vs Condon-Shortley %.2e  [%s]\n",
                    orb_label(oa) * orb_label(ob), F0, lb, Gl, worst, status)
            worst > tol && (ok = false)
        end
    end

    println(ok ? "===== Autoverificacion: TODO OK =====" :
                 "===== Autoverificacion: HAY FALLAS, no usar los numeros =====")
    return ok
end

# -----------------------------------------------------------------------------
#  Breit-Pauli (algebra identica a la de np2_toy_ci.jl, ya verificada vs limite jj)
# -----------------------------------------------------------------------------
function breit_pauli_levels(E3P, E1D, E1S, C3P, C1D, C1S, zeta)
    H0 = [E3P - zeta * C3P^2      sqrt(2.0) * zeta * C3P * C1S;
          sqrt(2.0) * zeta * C3P * C1S              E1S]
    e0 = eigvals(Symmetric(H0))

    E1 = E3P - 0.5 * zeta * C3P^2

    H2 = [E3P + 0.5 * zeta * C3P^2   -(1.0 / sqrt(2.0)) * zeta * C3P * C1D;
          -(1.0 / sqrt(2.0)) * zeta * C3P * C1D            E1D]
    f2 = eigen(Symmetric(H2))

    E3P0, E1S0 = e0[1], e0[2]
    E3P2, E1D2 = f2.values[1], f2.values[2]
    g_eff = f2.vectors[1, 1]^2 * 1.5 + f2.vectors[2, 1]^2 * 1.0

    return ([0.0,
             (E1 - E3P0) * AU2CM,
             (E3P2 - E3P0) * AU2CM,
             (E1D2 - E3P0) * AU2CM,
             (E1S0 - E3P0) * AU2CM], g_eff)
end

# -----------------------------------------------------------------------------
#  Driver
# -----------------------------------------------------------------------------
function run_full_ci(element_name::String, result_file::String;
                     n_per_l::Dict{Int,Int} = Dict(0 => 6, 1 => 6, 2 => 6, 3 => 4),
                     zeta_np::Union{Nothing,Float64} = nothing,
                     selftest::Bool = true, verbose::Bool = true,
                     ws = nothing)

    haskey(NP2_VALENCE_N, element_name) || error("Elemento no soportado: $element_name")
    n_val = NP2_VALENCE_N[element_name]

    data = load(result_file)
    R_max = data["R_max"]
    orbitals = data["orbitals"]
    alpha_d = get(data, "alpha_d", 0.0)
    r_c = get(data, "r_c", 1.0)

    if zeta_np === nothing
        key = "P_$(n_val)p"
        haskey(data, key) || error("$(result_file) no contiene $(key).")
        zeta_np = compute_zeta(data["R_grid"], data["V_eff"], data[key])
    end

    Z, N_el, K_ord, gamma = if element_name == "Carbon"
        (6.0, 100, 7, 2.5)
    elseif element_name == "Silicon"
        (14.0, 100, 7, 2.5)
    elseif element_name == "Germanium"
        (32.0, 300, 8, 3.0)
    else
        (50.0, 500, 8, 4.0)
    end
    # El workspace se puede compartir a lo largo de un estudio de convergencia: los pools
    # de tamanos crecientes estan ANIDADOS (Gram-Schmidt es secuencial y extract_virtuals
    # asigna el mismo pseudo_n independientemente de cuantos virtuales se pidan), y la
    # clave del cache es (n, l, k). Asi el espacio grande no recalcula los R^k del chico,
    # y length(ws.rk_cache) sigue siendo el numero de R^k distintos que ese espacio pide.
    if ws === nothing
        ws = K_ord == 7 ?
            cached_init_scf_workspace(R_max, N_el, Val(7), Z; γ = gamma, alpha_d = alpha_d, r_c = r_c) :
            cached_init_scf_workspace(R_max, N_el, Val(8), Z; γ = gamma, alpha_d = alpha_d, r_c = r_c)
    end

    val_orb = orbitals[end]
    core_orbs = orbitals[1:end-1]

    println("="^70)
    println(" CI de valencia completo -- $(element_name)")
    @printf(" archivo: %s | alpha_d = %.3f, r_c = %.3f\n", result_file, alpha_d, r_c)
    @printf(" zeta_%dp = %.8f Ha (%.3f cm^-1)\n", n_val, zeta_np, zeta_np * AU2CM)
    println(" espacio activo: " * join(["l=$l:$(n_per_l[l])" for l in sort(collect(keys(n_per_l)))], "  "))
    println("="^70)

    pool = build_orbital_pool(ws, core_orbs, val_orb, Z; n_per_l = n_per_l)
    val_idx = findfirst(o -> o.l == val_orb.l && o.n == val_orb.n, pool)
    val_idx === nothing && error("El orbital de valencia no quedo en el pool.")
    verbose && println("Pool (", length(pool), " orbitales): ",
                       join([orb_label(o) for o in pool], " "))

    if selftest && !engine_selftest(ws, pool, core_orbs, val_idx)
        error("Autoverificacion fallida: no se continua.")
    end

    results = Dict{String,Any}()
    for (L, S, name) in ((1, 1, "3P"), (2, 0, "1D"), (0, 0, "1S"))
        csfs = generate_pair_csfs(pool, L, S)
        idx = findfirst(c -> c.i == val_idx && c.j == val_idx, csfs)
        idx === nothing && error("El CSF de referencia np^2 no aparece en el bloque $name.")
        verbose && @printf("\nBloque %s: %d CSFs (referencia en la posicion %d)\n",
                           name, length(csfs), idx)
        H = build_pair_hamiltonian(ws, pool, core_orbs, csfs; verbose = verbose)
        F = eigen(Symmetric(H))
        E, C, root = pick_reference_root(F.values, F.vectors, idx)
        root != 1 && @warn "$name: la raiz de mayor peso np^2 es la #$root, no la mas baja."
        results[name] = (E = E, C = C, root = root, ref = H[idx, idx],
                         n = length(csfs), evals = F.values)
    end

    E3P, E1D, E1S = results["3P"].E, results["1D"].E, results["1S"].E
    C3P, C1D, C1S = results["3P"].C, results["1D"].C, results["1S"].C
    E_corr = E3P - results["3P"].ref

    levels, g_eff = breit_pauli_levels(E3P, E1D, E1S, C3P, C1D, C1S, zeta_np)

    println("\n--- Resultados CI ---")
    @printf("  E(3P) = %.10f Ha   (referencia %.10f, E_corr = %.6e Ha = %.1f cm^-1)\n",
            E3P, results["3P"].ref, E_corr, E_corr * AU2CM)
    @printf("  E(1D) = %.10f Ha\n", E1D)
    @printf("  E(1S) = %.10f Ha\n", E1S)
    @printf("  pesos np^2: C(3P) = %.6f  C(1D) = %.6f  C(1S) = %.6f\n", C3P, C1D, C1S)
    @printf("  1D - 3P = %10.2f cm^-1\n", (E1D - E3P) * AU2CM)
    @printf("  1S - 3P = %10.2f cm^-1\n", (E1S - E3P) * AU2CM)
    @printf("  g_eff(3P2) = %.6f\n", g_eff)

    println("\n--- Espectro de multipletes (cm^-1, desde 3P_0) ---")
    for (nm, v) in zip(["3P_0", "3P_1", "3P_2", "1D_2", "1S_0"], levels)
        @printf("  %-5s : %12.2f\n", nm, v)
    end
    println("-"^70)

    return (levels = levels, E_3P = E3P, E_1D = E1D, E_1S = E1S,
            E_corr = E_corr, zeta = zeta_np, C = (C3P, C1D, C1S),
            g_eff = g_eff, n_csf = (results["3P"].n, results["1D"].n, results["1S"].n),
            roots = (results["3P"].root, results["1D"].root, results["1S"].root),
            F2 = get_cached_Rk!(ws, pool[val_idx], pool[val_idx],
                                pool[val_idx], pool[val_idx], 2),
            n_orb = length(pool), n_per_l = n_per_l, n_rk = length(ws.rk_cache), ws = ws)
end

"""
    convergence_study(element, file; sizes, lmax)

Sube el espacio activo y reporta como converge E_corr. Es el numero que hay que
llevar a la tesis: no "el CI da X", sino "E_corr converge a X con este espacio".
"""
function convergence_study(element::String, result_file::String;
                           sizes = [3, 4, 6, 8], lmax::Int = 3)
    rows = []
    times = Float64[]
    ws = nothing
    for m in sizes
        npl = Dict(l => (l == lmax ? max(2, m - 2) : m) for l in 0:lmax)
        @printf("\n\n##### %s -- espacio activo m = %d, lmax = %d #####\n", element, m, lmax)
        t0 = time()
        r = run_full_ci(element, result_file; n_per_l = npl,
                        selftest = (m == first(sizes)), verbose = false, ws = ws)
        push!(times, time() - t0)
        push!(rows, r)
        ws = r.ws          # cache de R^k compartida con el siguiente tamano
    end

    println("\n\n===== Convergencia de E_corr -- $element (lmax = $lmax) =====")
    @printf("%4s %6s %8s %9s %14s %12s %10s %9s %s\n",
            "m", "orb", "CSF 3P", "R^k", "E_corr (Ha)", "E_corr(cm-1)", "3P_2", "t (s)", "raices")
    for (r, t) in zip(rows, times)
        @printf("%4d %6d %8d %9d %14.8f %12.1f %10.2f %9.1f  %d/%d/%d\n",
                r.n_per_l[1], r.n_orb, r.n_csf[1], r.n_rk, r.E_corr,
                r.E_corr * AU2CM, r.levels[3], t,
                r.roots[1], r.roots[2], r.roots[3])
    end
    # Las raices son el indice del autovector con mayor peso np^2 dentro de su bloque.
    # Que dejen de ser 1/1/1 no es un bug: significa que el espacio activo ya mete
    # estados s^2 o Rydberg por debajo de la referencia, y hay que decirlo en la tesis.
    any(r -> r.roots != (1, 1, 1), rows) &&
        println("\n  AVISO: en algun tamano la raiz de mayor peso np^2 dejo de ser la mas baja")
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    # Empezar chico: el carbono es el mas barato (100 elementos, K=7).
    convergence_study("Carbon", "carbon_rohf_results_3P_R30.0.jld2"; sizes = [3, 4, 6], lmax = 2)
end
