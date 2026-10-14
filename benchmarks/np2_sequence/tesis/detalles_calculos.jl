# =============================================================================
#  Calculos intermedios para el documento de detalles (docs/chapters/_borradores/detalles.tex)
#
#  Diagnostico de solo lectura: no toca ningun .jld2, ni resultados_ci.toml, ni docs/tablas/.
#  Reutiliza las funciones de la etapa 3 y escribe fragmentos tabular en
#  docs/chapters/_borradores/calculos/ con los pasos que el manuscrito solo cita como resultado:
#
#    extrapolaciones_datos.tex / _resultados.tex   seccion 6.2.1 (geometrica y de potencia)
#    ley_potencias_datos.tex / _ajuste.tex         seccion 6.3 (exponentes frente a Z)
#    mezcla_perturbativa.tex                       seccion 3.3 (serie perturbativa frente a la
#                                                  diagonalizacion exacta de cada bloque J)
#
#  Uso: julia --project=. benchmarks/np2_sequence/tesis/detalles_calculos.jl   (un minuto)
# =============================================================================

include(joinpath(@__DIR__, "etapa3_tablas.jl"))

const DIR_DETALLES = joinpath(REPO, "docs", "chapters", "_borradores", "calculos")

function escribir_fragmento(nombre, cuerpo, proc)
    mkpath(DIR_DETALLES)
    ruta = joinpath(DIR_DETALLES, nombre)
    open(ruta, "w") do io
        println(io, "% GENERADO por benchmarks/np2_sequence/tesis/detalles_calculos.jl; no editar a mano.")
        println(io, "% Procedencia:")
        foreach(p -> println(io, "%   ", p), proc)
        print(io, cuerpo)
    end
    println("  fragmento -> ", ruta)
end

pct(x, ref) = @sprintf("%+.1f", 100 * (x - ref) / ref)
TERM_LATEX = Dict(4 => raw"$^1D_2$", 5 => raw"$^1S_0$")

# -----------------------------------------------------------------------------
#  Extrapolacion de las curvas de convergencia de los singletes
# -----------------------------------------------------------------------------
function fragmentos_extrapolacion()
    io1, io2, proc = IOBuffer(), IOBuffer(), String[]
    println(io1, raw"\begin{tabular}{llrrrrrr}", "\n", raw"    \toprule")
    linea(io1, "Curva", raw"$m_1, m_2, m_3$", raw"$E(m_1)$", raw"$E(m_2)$", raw"$E(m_3)$",
          raw"$d_1$", raw"$d_2$", raw"$r = d_2/d_1$")
    println(io1, raw"    \midrule")
    println(io2, raw"\begin{tabular}{lrrrrrrr}", "\n", raw"    \toprule")
    linea(io2, "Curva", "NIST", raw"$E(m_3)$ (\%)", raw"$E_\infty^{\text{geo}}$", raw"(\%)",
          raw"$p$", raw"$E_\infty^{\text{pot}}$", raw"(\%)")
    println(io2, raw"    \midrule")
    for el in ELEMENTOS
        es = curva_convergencia(el)
        length(es) >= 3 || continue
        u = es[end-2:end]
        ms = [e["m"] for e in u]
        for t in (4, 5)
            ys = [e["niveles_cm"][t] for e in u]
            d1, d2 = ys[2] - ys[1], ys[3] - ys[2]
            g, p = extrap_geometrica(ys), extrap_potencia(ms, ys)
            nist = NIST_NIVELES[el][t]
            nombre = "$(INFO[el].etiqueta) $(TERM_LATEX[t])"
            linea(io1, nombre, join(ms, ", "), fmt(ys[1], 1), fmt(ys[2], 1), fmt(ys[3], 1),
                  fmt(d1, 1), fmt(d2, 1), fmt(d2 / d1, 2))
            linea(io2, nombre, fmt(nist, 1), pct(ys[3], nist),
                  g === nothing ? "--" : fmt(g.lim, 1), g === nothing ? "--" : pct(g.lim, nist),
                  p === nothing ? "--" : fmt(p.p, 1), p === nothing ? "--" : fmt(p.lim, 1),
                  p === nothing ? "--" : pct(p.lim, nist))
        end
        foreach(e -> push!(proc, proc_ci(e)), u)
    end
    for io in (io1, io2)
        println(io, raw"    \bottomrule", "\n", raw"\end{tabular}")
    end
    push!(proc, "Niveles en cm^-1 desde 3P_0; las mismas funciones extrap_geometrica y extrap_potencia de la etapa 3.")
    escribir_fragmento("extrapolaciones_datos.tex", String(take!(io1)), proc)
    escribir_fragmento("extrapolaciones_resultados.tex", String(take!(io2)), proc)
end

# -----------------------------------------------------------------------------
#  Ajuste de ley de potencias frente a Z
# -----------------------------------------------------------------------------
function fragmentos_potencias(HF)
    Zs = [INFO[el].Z for el in ELEMENTOS]
    series = [raw"$\zeta_{\text{calc}}$"     => [HF[el].zeta * HA2CM for el in ELEMENTOS],
              raw"$\zeta_{\text{desnudo}}$"  => [ALFA^2 / 2 * INFO[el].Z * HF[el].r3 * HA2CM for el in ELEMENTOS],
              raw"$\zeta_{\text{NIST}}$"     => [zeta_nist(el).zeta for el in ELEMENTOS],
              raw"$\langle r^{-3}\rangle$"   => [HF[el].r3 for el in ELEMENTOS],
              raw"$F^2(np,np)$"              => [HF[el].F2 * HA2CM for el in ELEMENTOS]]
    io1, io2 = IOBuffer(), IOBuffer()
    println(io1, raw"\begin{tabular}{lrrrrrrr}", "\n", raw"    \toprule")
    linea(io1, "Elemento", raw"$Z$", raw"$\ln Z$", [first(s) for s in series]...)
    println(io1, raw"    \midrule")
    for (i, el) in enumerate(ELEMENTOS)
        linea(io1, INFO[el].etiqueta, string(Int(Zs[i])), fmt(log(Zs[i]), 3),
              [fmt(last(s)[i], first(s) == raw"$\langle r^{-3}\rangle$" ? 3 : 1) for s in series]...)
    end
    println(io1, raw"    \bottomrule", "\n", raw"\end{tabular}")
    println(io2, raw"\begin{tabular}{lrrrrrr}", "\n", raw"    \toprule")
    linea(io2, "Serie", raw"$p$", raw"$A$", raw"$R^2$", raw"C\,$\to$\,Si", raw"Si\,$\to$\,Ge", raw"Ge\,$\to$\,Sn")
    println(io2, raw"    \midrule")
    for (nombre, ys) in series
        a = ajuste_potencia_Z(Zs, ys)
        loc = exponentes_locales(Zs, ys)
        linea(io2, nombre, fmt(a.p, 3), @sprintf("%.4g", a.A), fmt(a.R2, 4), [fmt(x, 2) for x in loc]...)
    end
    println(io2, raw"    \bottomrule", "\n", raw"\end{tabular}")
    proc = [proc_scf(archivo_scf(el, ESTADO)) for el in ELEMENTOS]
    push!(proc, "zeta y F^2 en cm^-1, <r^-3> en a0^-3; zeta_NIST de analisis_zeta.jl. Ajuste: ajuste_potencia_Z de la etapa 3.")
    escribir_fragmento("ley_potencias_datos.tex", String(take!(io1)), proc)
    escribir_fragmento("ley_potencias_ajuste.tex", String(take!(io2)), proc)
end

# -----------------------------------------------------------------------------
#  Serie perturbativa frente a diagonalizacion exacta en los bloques J = 0 y J = 2
# -----------------------------------------------------------------------------
"""
Bloque 2x2 [a x; x b] con b > a. Devuelve el cociente z = 2|x|/delta, el peso del segundo estado
en el autovector bajo a primer orden (x/delta)^2 y exacto, y el desplazamiento del nivel bajo a
segundo orden (-x^2/delta) y exacto.
"""
function bloque(a, b, x)
    δ = b - a
    z = 2 * abs(x) / δ
    w_exacto = (1 - 1 / sqrt(1 + z^2)) / 2
    return (δ = δ, x = x, z = z, w1 = (x / δ)^2, w = w_exacto,
            e2 = -x^2 / δ, e = -(sqrt(δ^2 / 4 + x^2) - δ / 2))
end

function fragmento_mezcla(HF)
    io, proc = IOBuffer(), String[]
    println(io, raw"\begin{tabular}{llrrrrrrr}", "\n", raw"    \toprule")
    linea(io, "Elemento", "Bloque", raw"$\delta$", raw"$x$", raw"$2|x|/\delta$",
          raw"$w^{(1)}$ (\%)", raw"$w$ (\%)", raw"$E^{(2)}$", raw"$E_{\text{exacto}}$")
    for el in ELEMENTOS
        e = ci_de(el, ESTADO, 0.0)
        e === nothing && continue
        ζ = HF[el].zeta * HA2CM
        E3P, E1D, E1S = e["E_3P_Ha"] * HA2CM, e["E_1D_Ha"] * HA2CM, e["E_1S_Ha"] * HA2CM
        C3P, C1D, C1S = e["pesos_np2"]
        # Los mismos elementos que breit_pauli_levels (np2_ci_full.jl), con los coeficientes CI.
        b2 = bloque(E3P + ζ * C3P^2 / 2, E1D, -ζ * C3P * C1D / sqrt(2))
        b0 = bloque(E3P - ζ * C3P^2, E1S, sqrt(2) * ζ * C3P * C1S)
        println(io, raw"    \midrule")
        for (k, (J, b)) in enumerate(((2, b2), (0, b0)))
            linea(io, k == 1 ? INFO[el].etiqueta : "", "\$J = $J\$", fmt(b.δ, 1), fmt(b.x, 2), fmt(b.z, 4),
                  fmt(100 * b.w1, 3), fmt(100 * b.w, 3), fmt(b.e2, 2), fmt(b.e, 2))
        end
        push!(proc, proc_ci(e))
    end
    println(io, raw"    \bottomrule", "\n", raw"\end{tabular}")
    push!(proc, "cm^-1. zeta calculado (3P, alpha_d = 0) y energias de los terminos del CI en el m de produccion. " *
                "delta = separacion entre los elementos diagonales; x = elemento fuera de la diagonal.")
    escribir_fragmento("mezcla_perturbativa.tex", String(take!(io)), proc)
end

function main_detalles()
    HF = Dict(el => propiedades_hf(el) for el in ELEMENTOS)
    fragmentos_extrapolacion()
    fragmentos_potencias(HF)
    fragmento_mezcla(HF)
    isempty(FALTAS) || println("Faltan: ", join(sort(collect(FALTAS)), "; "))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_detalles()
end
