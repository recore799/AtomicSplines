# =============================================================================
#  Diagnostico: energia ROHF con autovalores frente a cociente de Rayleigh
#
#  Los *_rohf.jl calculan E = sum_i (occ_i/2)(h_ii + eps_i) con eps_i el autovalor de la
#  Fock de su canal. El orbital de valencia, sin embargo, pasa despues por Gram-Schmidt
#  contra el core del mismo l, asi que deja de ser ese autovector y la expresion deja de
#  ser un funcional de los orbitales. Lo consistente es eps_i -> <c_i|F|c_i>.
#
#  Aqui NO se corre ningun SCF: el punto fijo lo definen las Fock, no la formula de la
#  energia, asi que basta reconstruir la Fock de cada canal con los orbitales guardados en
#  el .jld2 y evaluar las dos expresiones. Tres autoverificaciones:
#    1. la expresion con autovalores reproduce el E_total guardado;
#    2. los orbitales de core son autovectores de su canal: <c|F|c> = eps;
#    3. el autovalor de la Fock de valencia reconstruida reproduce el eps guardado.
#
#  Uso: julia --project=. benchmarks/np2_sequence/tesis/diagnostico_rayleigh.jl [--vpol]
#       Segundos por archivo. Con --vpol incluye tambien los barridos de V_pol.
# =============================================================================

isdefined(Main, :CABECERA_CI) || include(joinpath(@__DIR__, "comun.jl"))

"""
    focks_canal(ws, orbitals, estado) -> (F_core::Dict{Int,Matrix}, F_val)

Fock crudas (sin level shift ni extrapolacion) con la receta de los cuatro *_rohf.jl:
cada canal de core ve el core completo mas el promedio esferico de la capa abierta, y la
capa abierta ve el core mas su propio hueco de Fermi, con el coeficiente de intercambio
intra-capa del estado (2/25 promedio de configuracion, 5/25 termino 3P).
"""
function focks_canal(ws, orbitals, estado::String)
    n = ws.basis.num_splines
    core, val = orbitals[1:end-1], orbitals[end]
    coeff_k2 = estado == "av" ? 2.0 / 25.0 : 5.0 / 25.0
    H(l) = ws.T .+ ws.V .+ (l * (l + 1) / 2.0) .* ws.R_inv2

    build_total_J_matrix!(ws, core)
    J_core = copy(ws.J)
    J_esf = zeros(n, n)
    build_specific_J_matrix!(ws, J_esf, val, 2.0)

    F_core = Dict{Int,Matrix{Float64}}()
    for l in sort(unique(o.l for o in core))
        K_core = zeros(n, n); assemble_K_matrix!(ws, K_core, l, core)
        K_val  = zeros(n, n); assemble_K_matrix!(ws, K_val, l, [val])
        F_core[l] = H(l) .+ J_core .+ J_esf .- K_core .- K_val
    end

    K_core_p = zeros(n, n); assemble_K_matrix!(ws, K_core_p, 1, core)
    J_intra  = zeros(n, n); build_specific_J_matrix!(ws, J_intra, val, 1.0)
    K_intra  = zeros(n, n); build_specific_K_matrix!(ws, K_intra, val, 2, coeff_k2)
    F_val = H(1) .+ J_core .+ J_intra .- K_core_p .- K_intra
    return F_core, F_val
end

"""
    energias_rohf(ws, data, estado) -> NamedTuple

Energia con autovalores (recalculada), energia con cociente de Rayleigh y los tres
controles. `s2` es el peso que Gram-Schmidt le quita al autovector de valencia: la suma
de sus traslapes al cuadrado con el core del mismo l.
"""
function energias_rohf(ws, data, estado::String)
    orbitals = data["orbitals"]
    F_core, F_val = focks_canal(ws, orbitals, estado)
    val = orbitals[end]
    E_eps = 0.0
    E_ray = 0.0
    peor_core = 0.0
    for o in orbitals
        o.occ > 0 || continue
        c = o.coeffs
        h = dot(c, (ws.T .+ ws.V .+ (o.l * (o.l + 1) / 2.0) .* ws.R_inv2) * c)
        eps_R = dot(c, (o === val ? F_val : F_core[o.l]) * c)
        E_eps += (o.occ / 2) * (h + o.energy)
        E_ray += (o.occ / 2) * (h + eps_R)
        o === val || (peor_core = max(peor_core, abs(eps_R - o.energy)))
    end

    n = ws.basis.num_splines
    act = 3:(n - 1)
    ev, V = eigen(Symmetric(F_val[act, act]), Symmetric(ws.S[act, act]))
    j = argmax(abs.(V' * (ws.S[act, act] * val.coeffs[act])))
    v = zeros(n)
    v[act] = V[:, j]
    v ./= sqrt(dot(v, ws.S * v))
    s2 = sum((dot(o.coeffs, ws.S * v)^2 for o in orbitals[1:end-1] if o.l == 1); init = 0.0)

    return (E_guardada = data["E_total"], E_eps = E_eps, E_ray = E_ray,
            eps_val = val.energy, eps_reconstruido = ev[j], peor_core = peor_core, s2 = s2)
end

function main_rayleigh()
    man = manifiesto()
    con_vpol = "--vpol" in ARGS
    filas = []
    for (archivo, r) in sort(collect(man); by = p -> (orden_elemento(p[2]["elemento"]),
                                                      p[2]["alpha_d"], p[2]["estado"]))
        (r["alpha_d"] == 0.0 || con_vpol) || continue
        el = r["elemento"]
        verificar_malla(el, archivo, man)
        ws = workspace(el; alpha_d = r["alpha_d"])
        push!(filas, (archivo, energias_rohf(ws, load(ruta_np2(archivo)), r["estado"])))
    end

    println("\n", "="^118)
    println(" Energia ROHF: autovalores (lo que guardan los .jld2) frente a cociente de Rayleigh")
    println("="^118)
    @printf("%-46s %17s | %9s %9s %9s | %12s %9s\n", "archivo", "E guardada (Ha)",
            "ctrl 1", "ctrl 2", "ctrl 3", "E_ray-E_eps", "s2 (GS)")
    println("-"^118)
    for (archivo, e) in filas
        @printf("%-46s %17.8f | %9.1e %9.1e %9.1e | %12.3e %9.2e\n", archivo, e.E_guardada,
                e.E_eps - e.E_guardada, e.peor_core, e.eps_reconstruido - e.eps_val,
                e.E_ray - e.E_eps, e.s2)
    end
    println("-"^118)
    println(" ctrl 1: E con autovalores recalculada - E guardada.  ctrl 2: peor |<c|F|c> - eps| del core.")
    println(" ctrl 3: autovalor de la Fock de valencia reconstruida - eps guardado.")
    println(" Los tres tienen que ser chicos frente a E_ray - E_eps para que la comparacion valga.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_rayleigh()
end
