# =============================================================================
#  Diagnostico: casi-degeneracion ns^2 np^2 <-> np^4 con los orbitales HF guardados
#
#  np2_ci_full.jl correla solo la pareja np^2, con el ns^2 congelado. Este diagnostico mide
#  cuanto moveria los terminos la mezcla con np^4 en el modelo mas simple: dos configuraciones
#  por termino LS, con orbitales HF y core interno (todo menos ns y np) congelados,
#      H = [ 0  V ; V  Delta ],   V(3P) = V(1D) = G1(ns,np)/3,   V(1S) = 2 G1(ns,np)/3.
#  Los factores salen de <p^2 1S|1/r12|s^2 1S> = -G1/sqrt(3) por la amplitud de crear un par 1S
#  sobre p^2 (cuasiespin: 1/sqrt(3) con seniority 2 y 2/sqrt(3) con seniority 0). Delta =
#  E_av(np^4) - E_av(ns^2 np^2) es la misma para los tres terminos, porque np^4 y np^2 son capas
#  conjugadas con los mismos coeficientes de F^2: en este modelo el 1D - 3P no cambia.
#
#  Controles: las energias de orbital de ns y np guardadas en el .jld2 contra I + interacciones
#  promedio calculadas con las mismas integrales. Si pasan, las convenciones son las del SCF.
#
#  NO es una correccion para sumar al CI de pareja: los dos efectos no son aditivos
#  (docs/claude/REVISION-RESULTADOS-2026-09-13.md, seccion 7). Mide el tamano del efecto.
#
#  Uso: julia --project=. benchmarks/np2_sequence/tesis/diagnostico_casi_degeneracion.jl   (segundos)
# =============================================================================

isdefined(Main, :CABECERA_CI) || include(joinpath(@__DIR__, "comun.jl"))

"""
    casi_degeneracion(el, estado, alpha_d) -> NamedTuple

Integrales, controles y resultado del modelo de dos configuraciones sobre el .jld2 registrado.
Todo en Ha; `baja_S` es cuanto baja la separacion 1S - 3P.
"""
function casi_degeneracion(el, estado = ESTADO, alpha_d = 0.0)
    archivo = archivo_scf(el, estado, alpha_d)
    registrado(archivo, manifiesto()) || error("$archivo no esta registrado.")
    orbs = load(ruta_np2(archivo))["orbitals"]
    ws = workspace(el; alpha_d = alpha_d)
    p = orbs[end]
    s = only(filter(o -> o.l == 0 && o.n == INFO[el].n_val, orbs))
    interno = [o for o in orbs[1:end-1] if o !== s]
    nb = ws.basis.num_splines
    build_total_J_matrix!(ws, interno)
    J = copy(ws.J)
    function h(l)
        K = zeros(nb, nb)
        assemble_K_matrix!(ws, K, l, interno)
        return ws.T .+ ws.V .+ (l * (l + 1) / 2) .* ws.R_inv2 .+ J .- K
    end
    I_s = dot(s.coeffs, h(0) * s.coeffs)
    I_p = dot(p.coeffs, h(1) * p.coeffs)
    F0ss = compute_Rk(ws, s, s, s, s, 0)
    F0sp = compute_Rk(ws, s, p, s, p, 0)
    G1   = compute_Rk(ws, s, p, p, s, 1)
    F0pp = compute_Rk(ws, p, p, p, p, 0)
    F2pp = compute_Rk(ws, p, p, p, p, 2)
    k2 = estado == "av" ? 2 / 25 : 5 / 25
    ctrl_s = I_s + F0ss + 2F0sp - G1 / 3 - s.energy
    ctrl_p = I_p + 2F0sp - G1 / 3 + F0pp - k2 * F2pp - p.energy
    Δ = 2I_p - 2I_s + 5 * (F0pp - (2 / 25) * F2pp) - F0ss - 4 * (F0sp - G1 / 6)
    function dos_configuraciones(V)
        λ = (Δ - sqrt(Δ^2 + 4V^2)) / 2
        return (baja = -λ, peso = (λ / V)^2 / (1 + (λ / V)^2))
    end
    P, S = dos_configuraciones(G1 / 3), dos_configuraciones(2G1 / 3)
    return (archivo = archivo, G1 = G1, F2 = F2pp, Delta = Δ, ctrl_s = ctrl_s, ctrl_p = ctrl_p,
            peso_P = P.peso, peso_S = S.peso, baja_S = S.baja - P.baja)
end

function main_casi_degeneracion()
    ci = leer_ci()
    buscar(el, estado, a) = buscar_ci(ci, archivo_scf(el, estado, a), m_produccion(el))
    casos = [(el, 0.0) for el in ELEMENTOS]
    for el in ELEMENTOS
        haskey(BARRIDO_ALPHA, el) || continue
        sel = alpha_elegido(el, buscar)
        sel === nothing || push!(casos, (el, sel[1]))
    end
    sort!(casos; by = c -> (orden_elemento(c[1]), c[2]))

    println("\n", "="^134)
    println(" Casi-degeneracion ns^2 np^2 <-> np^4: dos configuraciones por termino, orbitales $(ESTADO) congelados")
    println("="^134)
    @printf("%-12s %8s %8s | %9s %9s | %8s %8s | %9s | %8s %8s -> %8s | %9s %9s\n",
            "caso", "G1", "Delta", "ctrl ns", "ctrl np", "np4 3P", "np4 1S", "baja S-P", "res 1D2", "res 1S0",
            "sumado", "HF+ND S-P", "HF D-P")
    println("-"^134)
    for (el, a) in casos
        r = casi_degeneracion(el, ESTADO, a)
        e = buscar(el, ESTADO, a)
        res4 = e === nothing ? NaN : e["niveles_cm"][4] - NIST_NIVELES[el][4]
        res5 = e === nothing ? NaN : e["niveles_cm"][5] - NIST_NIVELES[el][5]
        bS = r.baja_S * HA2CM
        @printf("%-12s %8.5f %8.4f | %9.1e %9.1e | %7.2f%% %7.2f%% | %9.0f | %+8.0f %+8.0f -> %+8.0f | %9.0f %9.0f\n",
                a > 0 ? @sprintf("%s ad=%.2f", el, a) : el, r.G1, r.Delta, r.ctrl_s, r.ctrl_p,
                100 * r.peso_P, 100 * r.peso_S, bS, res4, res5, res5 - bS, 0.6 * r.F2 * HA2CM - bS, 0.24 * r.F2 * HA2CM)
    end
    println("-"^134)
    println(" G1 y Delta en Ha; lo demas en cm^-1. res = CI de pareja - NIST en el m de produccion.")
    println(" 'sumado': residuo del 1S0 si se sumara la baja al CI de pareja. 'HF+ND S-P' = 0.6 F2 - baja, la")
    println(" separacion 1S - 3P con solo la casi-degeneracion; 'HF D-P' = 0.24 F2, que el modelo no cambia.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_casi_degeneracion()
end
