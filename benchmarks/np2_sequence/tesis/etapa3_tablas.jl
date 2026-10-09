# =============================================================================
#  Etapa 3: tablas del capitulo de resultados (segundos)
#
#  Lee los .jld2 registrados, resultados_ci.toml y las trazas de C-DIIS, y escribe en
#  docs/tablas/ un fragmento `tabular` por tabla, con su procedencia en comentarios, mas
#  valores_texto.md con los numeros que el texto cita fuera de las tablas. Los fragmentos
#  no llevan leyenda ni etiqueta; el capitulo los envuelve:
#      \begin{table} \centering \input{tablas/<nombre>} \caption{...} \label{...} \end{table}
#  Si falta un resultado (CI sin correr, SCF sin registrar) la celda sale como "--" y la
#  lista de faltantes lo dice al final: la etapa no inventa numeros.
#
#  Uso: julia --project=. benchmarks/np2_sequence/tesis/etapa3_tablas.jl
#       ... --ci /ruta/prueba.toml --m 3 --salida /ruta/dir      (para pruebas)
# =============================================================================

include(joinpath(@__DIR__, "comun.jl"))
isdefined(Main, :run_full_ci) || include(joinpath(NP2, "np2_ci_full.jl"))
include(joinpath(@__DIR__, "diagnostico_rayleigh.jl"))
include(joinpath(@__DIR__, "analisis_zeta.jl"))
include(joinpath(@__DIR__, "diagnostico_casi_degeneracion.jl"))
include(joinpath(@__DIR__, "diagnostico_estados_espurios.jl"))

const CI     = leer_ci(argumento("--ci", RESULTADOS_CI))
# --m fuerza un mismo tamano para todo (pruebas). Sin el, cada elemento usa su m de produccion y
# la sensibilidad se compara a M_SENSIBILIDAD (config.jl).
const M_ARG  = argumento("--m", nothing)
m_tabla(el)  = M_ARG === nothing ? m_produccion(el) : parse(Int, M_ARG)
m_sens()     = M_ARG === nothing ? M_SENSIBILIDAD : parse(Int, M_ARG)
const SALIDA = argumento("--salida", DIR_TABLAS)
const MAN    = manifiesto()
const FALTAS = Set{String}()

const TERMINOS = [raw"$^3P_0$", raw"$^3P_1$", raw"$^3P_2$", raw"$^1D_2$", raw"$^1S_0$"]
const TERM_TXT = ["3P_0", "3P_1", "3P_2", "1D_2", "1S_0"]

fmt(x, d) = x === nothing ? "--" : @sprintf("%.*f", d, x)
fmtS(x, d) = x === nothing ? "{--}" : @sprintf("%.*f", d, x)     # columnas S de siunitx
# Notacion cientifica sin ceros de relleno en el exponente (3.6e-8, no 3.6e-08), como la tesis.
cient(x) = replace(@sprintf("%.1E", x), r"E([+-])0*(\d)" => s"e\1\2")
# La misma notacion para columnas que no son de siunitx: $4.0\times10^{-12}$.
cient_tex(x) = replace(cient(x), r"e([+-]?)(\d+)" => s"\\times10^{\1\2}") |> t -> "\$" * replace(t, "{+" => "{") * "\$"
linea(io, partes...) = println(io, "    ", join(partes, " & "), " \\\\")
etiqueta_Z(el) = "$(INFO[el].etiqueta) (\$Z=$(Int(INFO[el].Z))\$)"
proc_scf(archivo) = "$archivo (sha256 $(corto(sha256_de(ruta_np2(archivo)))))"
proc_ci(e) = "CI $(e["archivo"]) m=$(e["m"]) (sha256 $(corto(e["sha256_entrada"])), commit $(e["commit"]))"

function ci_de(el, estado, a; m::Int = estado == ESTADO ? m_tabla(el) : m_sens())
    archivo = archivo_scf(el, estado, a)
    e = registrado(archivo, MAN) ? buscar_ci(CI, archivo, m) : nothing
    e === nothing && push!(FALTAS, "CI de $archivo con m = $m")
    return e
end

# -----------------------------------------------------------------------------
#  Propiedades Hartree-Fock (sin CI) de cada elemento, sobre el .jld2 de ESTADO
# -----------------------------------------------------------------------------
function propiedades_hf(el, estado = ESTADO, a = 0.0)
    archivo = archivo_scf(el, estado, a)
    registrado(archivo, MAN) || error("$archivo no esta registrado: correr la etapa 1.")
    verificar_malla(el, archivo, MAN)
    data = load(ruta_np2(archivo))
    ws = workspace(el; alpha_d = a)
    orbs = data["orbitals"]
    v = orbs[end]
    T = sum(o.occ * dot(o.coeffs, (ws.T .+ (o.l * (o.l + 1) / 2.0) .* ws.R_inv2) * o.coeffs)
            for o in orbs if o.occ > 0)
    E = ENERGIA_HF == :rayleigh ? energias_rohf(ws, data, estado).E_ray : data["E_total"]
    R, Veff = data["R_grid"], data["V_eff"]
    Vrad = Veff .+ 1.0 ./ R .^ 2          # l = 1: l(l+1)/2r^2 = 1/r^2
    i = argmin(Vrad)
    return (archivo = archivo, E = E, T = T, virial = -(E - T) / T, eps = v.energy,
            r3 = dot(v.coeffs, ws.R_inv3 * v.coeffs),
            F0 = compute_Rk(ws, v, v, v, v, 0), F2 = compute_Rk(ws, v, v, v, v, 2),
            zeta = compute_zeta(R, Veff, data["P_$(INFO[el].n_val)p"]),
            vrad_min = Vrad[i], r_min = R[i])
end


alpha_elegido(el) = alpha_elegido(el, ci_de)
curva_convergencia(el) = curva_convergencia(CI, MAN, el)

# -----------------------------------------------------------------------------
#  Tablas
# -----------------------------------------------------------------------------
function tabla_energia(HF)
    io = IOBuffer(); proc = String[]
    println(io, raw"\begin{tabular}{l l S[table-format=-4.8] S[table-format=-4.10] S[table-format=1.1e-1]}")
    println(io, raw"    \toprule")
    linea(io, raw"\textbf{Elemento}", raw"\textbf{Magnitud}", raw"{\textbf{Este trabajo}}",
          raw"{\textbf{Froese Fischer}}", raw"{$|\Delta|$}")
    for el in ELEMENTOS
        h, ref = HF[el], FF[el]
        println(io, raw"    \midrule")
        println(io, "    \\multirow{3}{*}{$(INFO[el].etiqueta)}")
        for (nombre, x, r) in ((raw"$E$ ($E_h$)", h.E, ref.E), (raw"$T$ ($E_h$)", h.T, ref.T),
                               (raw"$-V/T$", h.virial, ref.virial))
            d = r === nothing ? "{--}" : cient(abs(x - parse(Float64, r)))
            linea(io, "", nombre, @sprintf("%.8f", x), something(r, "{--}"), d)
        end
        push!(proc, proc_scf(h.archivo) * (ENERGIA_HF == :rayleigh ? ", E por cociente de Rayleigh" : ""))
    end
    println(io, raw"    \bottomrule", "\n", raw"\end{tabular}")
    push!(proc, "Referencias: Froese Fischer (1977), copiadas en tesis/config.jl.")
    return String(take!(io)), proc
end

"""F^0, F^2 y <r^-3> del orbital de valencia frente a Froese Fischer, en un solo cuadro."""
function tabla_radiales(HF)
    io = IOBuffer(); proc = String[]
    println(io, raw"\begin{tabular}{l l S[table-format=1.8] S[table-format=1.8] S[table-format=1.1e-1]}")
    println(io, raw"    \toprule")
    linea(io, raw"\textbf{Elemento}", raw"\textbf{Magnitud}", raw"{\textbf{Este trabajo}}",
          raw"{\textbf{Froese Fischer}}", raw"{$|\Delta|$}")
    for el in ELEMENTOS
        h, n, ref = HF[el], INFO[el].n_val, FF[el]
        println(io, raw"    \midrule")
        println(io, "    \\multirow{3}{*}{$(INFO[el].etiqueta)}")
        for (nombre, x, r, dec) in (("\$F^0($(n)p,$(n)p)\$ (\$E_h\$)", h.F0, ref.F0, 8),
                                    ("\$F^2($(n)p,$(n)p)\$ (\$E_h\$)", h.F2, ref.F2, 8),
                                    ("\$\\langle r^{-3}\\rangle_{$(n)p}\$ (\$a_0^{-3}\$)", h.r3, ref.r3, 6))
            d = r === nothing ? "{--}" : cient(abs(x - parse(Float64, r)))
            linea(io, "", nombre, @sprintf("%.*f", dec, x), something(r, "{--}"), d)
        end
        push!(proc, proc_scf(h.archivo))
    end
    println(io, raw"    \bottomrule", "\n", raw"\end{tabular}")
    push!(proc, "Referencias: Froese Fischer (1977), copiadas en tesis/config.jl.")
    return String(take!(io)), proc
end

# Error relativo con signo frente a una referencia, en por ciento.
err_rel(v, ref) = (v === nothing || ref == 0.0) ? "--" : @sprintf("\$%+.1f\$\\%%", 100 * (v - ref) / ref)

"""Energía de ionización de los cuatro elementos: Koopmans, + CI de valencia, + V_pol, NIST."""
function tabla_ionizacion(HF)
    io = IOBuffer(); proc = String[]
    celda(x, nist) = x === nothing ? "--" : "$(fmt(x, 2)) ($(err_rel(x, nist)))"
    println(io, raw"\begin{tabular}{lcccc}")
    println(io, raw"    \toprule")
    linea(io, raw"\textbf{Elemento}", raw"\textbf{Koopmans}", raw"\textbf{HF+CI}",
          raw"\textbf{CI+$V_{\text{pol}}$}", raw"\textbf{NIST}")
    println(io, raw"    \midrule")
    for el in ELEMENTOS
        h, e0, nist = HF[el], ci_de(el, ESTADO, 0.0), NIST_IONIZACION[el]
        ip_hf = -h.eps * HA2EV
        ip_ci = e0 === nothing ? nothing : (-h.eps - e0["E_corr_Ha"]) * HA2EV
        push!(proc, proc_scf(h.archivo))
        e0 === nothing || push!(proc, proc_ci(e0))
        ip_vp = nothing
        if haskey(BARRIDO_ALPHA, el)
            sel = alpha_elegido(el)
            if sel !== nothing
                a, ev = sel
                eps_a = load(ruta_np2(ev["archivo"]))["orbitals"][end].energy
                ip_vp = (-eps_a - ev["E_corr_Ha"]) * HA2EV
                push!(proc, proc_scf(ev["archivo"]), proc_ci(ev) * @sprintf(" <- alpha_d = %.2f", a))
            end
        end
        linea(io, INFO[el].etiqueta, celda(ip_hf, nist), celda(ip_ci, nist), celda(ip_vp, nist), fmt(nist, 2))
    end
    println(io, raw"    \bottomrule", "\n", raw"\end{tabular}")
    push!(proc, "En eV, con el error relativo frente al NIST entre parentesis. HF+CI = -eps_np - E_corr: " *
                "el ion np^1 no tiene correlacion de pareja de valencia.")
    return String(take!(io)), proc
end

"""
Niveles de los elementos `els` frente al NIST, un bloque de filas por elemento. El 3P_0 es el
origen de la escala y no se imprime. Con `con_vpol`, se añaden las columnas CI+V_pol del alpha_d
elegido.
"""
function tabla_niveles(els, con_vpol::Bool)
    io = IOBuffer(); proc = String[]
    println(io, con_vpol ? raw"\begin{tabular}{llrrrrr}" : raw"\begin{tabular}{llrrr}")
    println(io, raw"    \toprule")
    enc = con_vpol ? [raw"\textbf{HF+CI}", raw"$\Delta$", raw"\textbf{CI+$V_{\text{pol}}$}", raw"$\Delta$", raw"\textbf{NIST}"] :
                     [raw"\textbf{HF+CI}", raw"$\Delta$", raw"\textbf{NIST}"]
    linea(io, raw"\textbf{Elemento}", raw"\textbf{Nivel}", enc...)
    for el in els
        e0 = ci_de(el, ESTADO, 0.0)
        e0 === nothing || push!(proc, proc_ci(e0))
        ev = nothing
        if con_vpol
            sel = alpha_elegido(el)
            if sel !== nothing
                ev = sel[2]
                push!(proc, proc_ci(ev) * @sprintf(" <- alpha_d elegido por %s", CRITERIO_ALPHA))
            end
        end
        println(io, raw"    \midrule")
        println(io, "    \\multirow{4}{*}{$(INFO[el].etiqueta)}")
        for t in 2:5
            ref = NIST_NIVELES[el][t]
            v0 = e0 === nothing ? nothing : e0["niveles_cm"][t]
            celdas = [fmt(v0, 1), err_rel(v0, ref)]
            if con_vpol
                vv = ev === nothing ? nothing : ev["niveles_cm"][t]
                append!(celdas, [fmt(vv, 1), err_rel(vv, ref)])
            end
            linea(io, "", TERMINOS[t], celdas..., fmt(ref, 1))
        end
    end
    println(io, raw"    \bottomrule", "\n", raw"\end{tabular}")
    push!(proc, "Niveles en cm^-1 desde el 3P_0 (origen de la escala, no se imprime). " *
                "Delta = error relativo frente al NIST, en por ciento.")
    return String(take!(io)), proc
end

"""Barrido de alpha_d de todos los elementos que lo tienen, con zeta en cm^-1."""
function tabla_barrido()
    io = IOBuffer(); proc = String[]
    println(io, raw"\begin{tabular}{lcrrrr}")
    println(io, raw"    \toprule")
    linea(io, raw"\textbf{Elemento}", raw"$\alpha_d$", raw"$\zeta_{np}$", raw"$^3P_1$", raw"$^3P_2$",
          CRITERIO_ALPHA == :solo_3P2 ? raw"Error $^3P_2$" : raw"Error rms")
    for el in ELEMENTOS
        haskey(BARRIDO_ALPHA, el) || continue
        sel = alpha_elegido(el)
        as = vcat(0.0, BARRIDO_ALPHA[el])
        println(io, raw"    \midrule")
        println(io, "    \\multirow{$(length(as))}{*}{$(INFO[el].etiqueta)}")
        for a in as
            e = ci_de(el, ESTADO, a)
            marca = (sel !== nothing && sel[1] == a) ? raw"$^{*}$" : ""
            if e === nothing
                linea(io, "", fmt(a, 2) * marca, "--", "--", "--", "--")
            else
                n = e["niveles_cm"]
                linea(io, "", fmt(a, 2) * marca, fmt(e["zeta_Ha"] * HA2CM, 1), fmt(n[2], 1), fmt(n[3], 1),
                      @sprintf("%.1f\\%%", 100 * costo_alpha(n, el)))
                push!(proc, proc_ci(e))
            end
        end
    end
    println(io, raw"    \bottomrule", "\n", raw"\end{tabular}")
    push!(proc, "zeta y niveles en cm^-1. El asterisco marca el alpha_d elegido por el criterio " *
                "$(CRITERIO_ALPHA); r_c = $(R_C) a0.")
    return String(take!(io)), proc
end

function tabla_convergencia()
    io = IOBuffer(); io2 = IOBuffer(); proc = String[]
    por_el = Dict(el => Dict(e["m"] => e for e in CI if e["archivo"] == archivo_scf(el, ESTADO, 0.0) &&
                                                      e["sha256_entrada"] == sha256_de(ruta_np2(archivo_scf(el, ESTADO, 0.0))))
                  for el in ELEMENTOS)
    ms = sort(unique(vcat([collect(keys(d)) for d in values(por_el)]...)))
    nel = length(ELEMENTOS)
    println(io, raw"\begin{tabular}{rr", "r"^nel, "}")
    println(io, raw"    \toprule")
    linea(io, "", "", "\\multicolumn{$nel}{c}{\$E_{\\text{corr}}\$ (mHa)}")
    println(io, "    \\cmidrule(lr){3-$(2 + nel)}")
    linea(io, raw"$m$", raw"CSF $^3P$", [el for el in ELEMENTOS]...)
    println(io, raw"    \midrule")
    println(io2, raw"\begin{tabular}{r", "r"^(2 * nel), "}")
    println(io2, raw"    \toprule")
    linea(io2, "", "\\multicolumn{$nel}{c}{Integrales \$R^k\$ en caché}", "\\multicolumn{$nel}{c}{Tiempo de pared (s)}")
    println(io2, "    \\cmidrule(lr){2-$(1 + nel)} \\cmidrule(lr){$(2 + nel)-$(1 + 2nel)}")
    linea(io2, raw"$m$", [el for el in ELEMENTOS]..., [el for el in ELEMENTOS]...)
    println(io2, raw"    \midrule")
    for m in ms
        ref = first(d[m] for d in values(por_el) if haskey(d, m))
        linea(io, string(m), string(ref["n_csf"][1]),
              [haskey(por_el[el], m) ? "\$$(fmt(1e3 * por_el[el][m]["E_corr_Ha"], 3))\$" : "--" for el in ELEMENTOS]...)
        linea(io2, string(m), [haskey(por_el[el], m) ? string(por_el[el][m]["n_rk"]) : "--" for el in ELEMENTOS]...,
              [haskey(por_el[el], m) ? fmt(por_el[el][m]["tiempo_s"], 0) : "--" for el in ELEMENTOS]...)
    end
    for io_ in (io, io2)
        println(io_, raw"    \bottomrule", "\n", raw"\end{tabular}")
    end
    for el in ELEMENTOS, (m, e) in sort(collect(por_el[el]); by = first)
        push!(proc, proc_ci(e))
    end
    push!(proc, "lmax = $(LMAX_CI). t es incremental (cada tamano reutiliza las R^k del anterior) salvo " *
                "donde resultados_ci.toml marca cache_heredada = false: primer tamano o corrida reanudada.")
    return String(take!(io)), String(take!(io2)), proc
end

function traza_diis(el, tag)
    ruta = joinpath(NP2, "diis_trace_$(INFO[el].prefijo)_$(tag).csv")
    isfile(ruta) || return nothing
    # Las trazas sin C-DIIS de Ge y Sn llevan comentarios de procedencia en la cabecera.
    filas = [split(l, ",") for l in readlines(ruta) if !isempty(strip(l)) && !startswith(l, "#") && !startswith(l, "iter")]
    E = strip(filas[end][2])
    decimales = occursin('.', E) ? length(split(E, '.')[2]) : 0
    return (iters = length(filas), E = parse(Float64, E), decimales = decimales, archivo = basename(ruta))
end

function tabla_cdiis()
    io = IOBuffer(); proc = String[]
    println(io, raw"\begin{tabular}{lrrrc}")
    println(io, raw"    \toprule")
    linea(io, "", raw"\multicolumn{2}{c}{\textbf{Iteraciones}}", "", "")
    println(io, raw"    \cmidrule(lr){2-3}")
    linea(io, raw"\textbf{Elemento}", raw"sin C-DIIS", raw"con C-DIIS", raw"\textbf{Factor}", raw"$|\Delta E|$ (Ha)")
    println(io, raw"    \midrule")
    for el in ELEMENTOS
        s, c = traza_diis(el, "sin_diis"), traza_diis(el, "con_diis")
        if s === nothing || c === nothing
            push!(FALTAS, "trazas de C-DIIS de $(INFO[el].prefijo)")
            linea(io, INFO[el].etiqueta, "--", "--", "--", "--")
            continue
        end
        # Una diferencia por debajo de la ultima cifra impresa en la traza no se puede afirmar.
        d = min(s.decimales, c.decimales)
        dE = abs(s.E - c.E)
        linea(io, INFO[el].etiqueta, string(s.iters), string(c.iters), @sprintf("%.1f", s.iters / c.iters),
              dE < 10.0^(-d) ? "\$<10^{-$d}\$" : cient_tex(dE))
        push!(proc, "$(s.archivo), $(c.archivo)")
    end
    println(io, raw"    \bottomrule", "\n", raw"\end{tabular}")
    push!(proc, "Trazas de diis_benchmark.jl (estado 3P, |dE| < 1e-10 Ha, mismo punto de partida). " *
                "Las sin C-DIIS de Ge y Sn se reconstruyeron del log impreso, con E a 8 decimales.")
    return String(take!(io)), proc
end

"""
Residuales de las trazas con C-DIIS: los de la ultima iteracion (conmutador crudo y proyectado)
y la primera iteracion en que el proyectado queda a menos del 10 % de su valor final, es decir,
donde toca su piso. `nothing` si la traza no existe.
"""
function residuales_diis(el)
    ruta = joinpath(NP2, "diis_trace_$(INFO[el].prefijo)_con_diis.csv")
    isfile(ruta) || return nothing
    filas = [parse.(Float64, split(l, ",")) for l in readlines(ruta)
             if !isempty(strip(l)) && !startswith(l, "#") && !startswith(l, "iter")]
    proy, crudo = filas[end][4], filas[end][5]
    piso = findfirst(f -> f[4] <= 1.1 * proy, filas)
    return (iters = length(filas), proy = proy, crudo = crudo, iter_piso = Int(filas[piso][1]))
end

"""
Base de B-splines de cada elemento: la malla de INFO (comun.jl), que `verificar_malla` contrasta
con RESULTS.toml, y el numero de splines leido de los coeficientes guardados en el .jld2.
"""
function tabla_parametros_base()
    io = IOBuffer(); proc = String[]
    println(io, raw"\begin{tabular}{lccccc}")
    println(io, raw"    \toprule")
    linea(io, raw"\textbf{Elemento}", raw"$R_{\text{max}}$ ($a_0$)", raw"\textbf{Intervalos}",
          raw"\textbf{Orden} $k$", raw"$\gamma$", raw"\textbf{Splines}")
    println(io, raw"    \midrule")
    for el in ELEMENTOS
        i = INFO[el]
        archivo = archivo_scf(el, ESTADO)
        verificar_malla(el, archivo, MAN)
        n = length(load(ruta_np2(archivo))["orbitals"][1].coeffs)
        linea(io, i.etiqueta, fmt(R_MAX, 0), string(i.N), string(i.K), fmt(i.gamma, 1), string(n))
        push!(proc, proc_scf(archivo))
    end
    println(io, raw"    \bottomrule", "\n", raw"\end{tabular}")
    push!(proc, "Malla de INFO (comun.jl), la de cada *_rohf.jl y de RESULTS.toml. Splines = longitud de los " *
                "coeficientes guardados; k es el orden (grado k - 1).")
    return String(take!(io)), proc
end

function tabla_zeta(HF, ZETA_SENS)
    io = IOBuffer(); proc = String[]
    println(io, raw"\begin{tabular}{lccccc}")
    println(io, raw"    \toprule")
    linea(io, raw"\textbf{Elemento}", raw"$(\alpha Z)^2$", raw"$\zeta_{\text{desnudo}}$",
          raw"$\zeta_{\text{calc}}$", raw"$\zeta_{\text{NIST}}$", raw"$\zeta_{\text{calc}}/\zeta_{\text{NIST}}$")
    println(io, raw"    \midrule")
    for el in ELEMENTOS
        h, Z = HF[el], INFO[el].Z
        zb = ALFA^2 / 2 * Z * h.r3 * HA2CM
        zc = h.zeta * HA2CM
        zn = zeta_nist(el).zeta
        linea(io, INFO[el].etiqueta, fmt((ALFA * Z)^2, 4), fmt(zb, 1), fmt(zc, 1), fmt(zn, 1), fmt(zc / zn, 3))
        push!(proc, proc_scf(h.archivo))
    end
    println(io, raw"    \bottomrule", "\n", raw"\end{tabular}")
    push!(proc, "zeta en cm^-1. zeta_desnudo = (alpha^2/2) Z <r^-3>. zeta_NIST: la zeta unica que reproduce " *
                "3P_2, 1D_2 y 1S_0 del NIST en Breit-Pauli de p^2 sin CI, con las energias LS libres.")
    return String(take!(io)), proc
end

"""
Escalas del regimen de acoplamiento: zeta frente a la separacion electrostatica
E(1D) - E(3P) = (6/25) F^2 de Hartree-Fock, y la razon de intervalos R del modelo y del NIST.
"""
function tabla_acoplamiento(HF)
    io = IOBuffer(); proc = String[]
    println(io, raw"\begin{tabular}{lccccc}")
    println(io, raw"    \toprule")
    linea(io, raw"\textbf{Elemento}", raw"$\zeta_{np}$", raw"$\frac{6}{25}F^2$",
          raw"$\zeta_{np}\big/\frac{6}{25}F^2$", raw"$R_{\text{modelo}}$", raw"$R_{\text{NIST}}$")
    println(io, raw"    \midrule")
    for el in ELEMENTOS
        h = HF[el]
        zc, dls = h.zeta * HA2CM, 6 / 25 * h.F2 * HA2CM
        n = NIST_NIVELES[el]
        e0 = ci_de(el, ESTADO, 0.0)
        Rmod = e0 === nothing ? nothing : (e0["niveles_cm"][3] - e0["niveles_cm"][2]) / e0["niveles_cm"][2]
        linea(io, INFO[el].etiqueta, fmt(zc, 1), fmt(dls, 0), fmt(zc / dls, 4), fmt(Rmod, 3), fmt((n[3] - n[2]) / n[2], 3))
        push!(proc, proc_scf(h.archivo))
        e0 === nothing || push!(proc, proc_ci(e0))
    end
    println(io, raw"    \bottomrule", "\n", raw"\end{tabular}")
    push!(proc, "zeta y (6/25)F^2 = E(1D) - E(3P) de Hartree-Fock en cm^-1. R = (3P_2 - 3P_1)/(3P_1 - 3P_0); " *
                "una zeta de un cuerpo en LS puro da R = 2.")
    return String(take!(io)), proc
end

# g de Lande de 3P_2 y 1D_2 puros con el g_s del electron libre. El NIST mide con el g_s real: con
# g_s = 2 el 3P_2 puro valdria 1.5 en vez de 1.50116, un corrimiento mayor que la mezcla del Ge.
const G_3P2 = 1.0 + (G_S - 1.0) / 2      # 1 + (g_s - 1)[J(J+1) - L(L+1) + S(S+1)]/[2J(J+1)]
const G_1D2 = 1.0
peso_1D(g) = (G_3P2 - g) / (G_3P2 - G_1D2)

"""g del 3P_2 con g_s real, a partir del g_eff de run_full_ci, que usa g_s = 2 (1.5 y 1.0)."""
function g_real(e)
    c2 = 2 * (1.5 - e["g_eff_3P2"])      # peso de 1D_2 en el nivel
    return G_3P2 * (1 - c2) + G_1D2 * c2
end

function tabla_lande()
    io = IOBuffer(); proc = String[]
    pct(x) = x === nothing ? "--" : @sprintf("%.2f", 100 * x)
    println(io, raw"\begin{tabular}{lcccccc}")
    println(io, raw"    \toprule")
    linea(io, "", raw"\multicolumn{3}{c}{\textbf{Factor $g$ del $^3P_2$}}", raw"\multicolumn{3}{c}{\textbf{Peso de $^1D_2$ (\%)}}")
    println(io, raw"    \cmidrule(lr){2-4} \cmidrule(lr){5-7}")
    linea(io, raw"\textbf{Elemento}", "HF+CI", raw"CI+$V_{\text{pol}}$", "NIST", "HF+CI", raw"CI+$V_{\text{pol}}$", "NIST")
    println(io, raw"    \midrule")
    for el in ELEMENTOS
        e0 = ci_de(el, ESTADO, 0.0)
        sel = haskey(BARRIDO_ALPHA, el) ? alpha_elegido(el) : nothing
        g0 = e0 === nothing ? nothing : g_real(e0)
        gv = sel === nothing ? nothing : g_real(sel[2])
        gn = NIST_LANDE_3P2[el]
        linea(io, INFO[el].etiqueta, fmt(g0, 5), fmt(gv, 5), gn === nothing ? "--" : string(gn),
              pct(g0 === nothing ? nothing : peso_1D(g0)), pct(gv === nothing ? nothing : peso_1D(gv)),
              pct(gn === nothing ? nothing : peso_1D(gn)))
        e0 === nothing || push!(proc, proc_ci(e0))
        sel === nothing || push!(proc, proc_ci(sel[2]))
    end
    println(io, raw"    \bottomrule", "\n", raw"\end{tabular}")
    push!(proc, "g del 3P_2 con g_s = $(G_S), sin correcciones relativistas ni diamagneticas (orden alpha^2). " *
                "Peso de 1D_2 en el nivel = (g(3P_2 puro) - g)/(g(3P_2 puro) - g(1D_2 puro)). NIST ASD ver. 5.12.")
    return String(take!(io)), proc
end

"""Ley de potencias y = A Z^p por minimos cuadrados sobre log-log, con su R^2."""
function ajuste_potencia_Z(Zs, ys)
    x, y = log.(Zs), log.(ys)
    mx, my = sum(x) / length(x), sum(y) / length(y)
    p = sum((x .- mx) .* (y .- my)) / sum((x .- mx) .^ 2)
    res = y .- (my .+ p .* (x .- mx))
    tot = y .- my
    return (p = p, A = exp(my - p * mx), R2 = 1 - sum(res .^ 2) / sum(tot .^ 2))
end

"""Exponente local entre elementos consecutivos: log(y2/y1) / log(Z2/Z1)."""
exponentes_locales(Zs, ys) = [log(ys[i + 1] / ys[i]) / log(Zs[i + 1] / Zs[i]) for i in 1:length(Zs) - 1]

"""Limite si los incrementos siguen en razon constante; nothing si la razon no esta en (0, 1)."""
function extrap_geometrica(y)
    d1, d2 = y[2] - y[1], y[3] - y[2]
    r = d2 / d1
    (0 < r < 1) || return nothing
    return (lim = y[3] + d2 * r / (1 - r), razon = r)
end

"""Ajuste exacto de E(m) = E_inf + A m^-p a tres puntos; nothing si no hay un p > 0 que lo cumpla."""
function extrap_potencia(m, y)
    q = (y[2] - y[1]) / (y[3] - y[2])
    f(p) = (m[1]^-p - m[2]^-p) / (m[2]^-p - m[3]^-p) - q
    lo, hi = 0.05, 20.0
    f(lo) * f(hi) < 0 || return nothing
    for _ in 1:200
        mid = (lo + hi) / 2
        f(lo) * f(mid) <= 0 ? (hi = mid) : (lo = mid)
    end
    p = (lo + hi) / 2
    A = (y[3] - y[2]) / (m[3]^-p - m[2]^-p)
    return (lim = y[3] - A * m[3]^-p, p = p)
end

function tabla_convergencia_singletes()
    io = IOBuffer(); proc = String[]
    curvas = Dict(el => Dict(e["m"] => e for e in curva_convergencia(el)) for el in ELEMENTOS)
    ms = sort(unique(vcat([collect(keys(c)) for c in values(curvas)]...)))
    println(io, raw"\begin{tabular}{r", "rr"^length(ELEMENTOS), "}")
    println(io, raw"    \toprule")
    linea(io, "", ["\\multicolumn{2}{c}{\\textbf{$(INFO[el].etiqueta)}}" for el in ELEMENTOS]...)
    linea(io, raw"$m$", repeat([raw"$^1D_2$", raw"$^1S_0$"], length(ELEMENTOS))...)
    println(io, raw"    \midrule")
    for m in ms
        linea(io, string(m), vcat([[haskey(curvas[el], m) ? fmt(curvas[el][m]["niveles_cm"][t], 1) : "--"
                                    for t in (4, 5)] for el in ELEMENTOS]...)...)
    end
    println(io, raw"    \midrule")
    linea(io, "NIST", vcat([[fmt(NIST_NIVELES[el][t], 1) for t in (4, 5)] for el in ELEMENTOS]...)...)
    println(io, raw"    \bottomrule", "\n", raw"\end{tabular}")
    for el in ELEMENTOS, m in sort(collect(keys(curvas[el])))
        push!(proc, proc_ci(curvas[el][m]))
    end
    push!(proc, "Niveles en cm^-1 desde 3P_0, lmax = $(LMAX_CI). Las extrapolaciones estan en valores_texto.md.")
    return String(take!(io)), proc
end

function valores_texto(HF, ZETA_SENS)
    io = IOBuffer()
    println(io, "# Valores citados en el texto\n")
    println(io, "GENERADO por `benchmarks/np2_sequence/tesis/etapa3_tablas.jl`; no editar a mano.")
    println(io, "Convencion de orbitales: $(ESTADO) (sensibilidad: $(ESTADO_SENSIBILIDAD) a m = $(m_sens())); ",
            "lmax = $(LMAX_CI); m de produccion: ", join(["$(el) $(m_tabla(el))" for el in ELEMENTOS], ", "), ".\n")
    println(io, "## Pozo del potencial radial efectivo V_rad = V_eff + 1/r^2\n")
    for el in ELEMENTOS
        @printf(io, "- %s: minimo %.2f Ha en r = %.3f a0\n", INFO[el].etiqueta, HF[el].vrad_min, HF[el].r_min)
    end
    println(io, "\n## zeta y F^2 (cm^-1)\n")
    for el in ELEMENTOS
        h = HF[el]
        @printf(io, "- %s: zeta(%s) = %.3f, zeta(%s) = %.3f (razon %.4f), F^2 = %.1f\n", INFO[el].etiqueta,
                ESTADO, h.zeta * HA2CM, ESTADO_SENSIBILIDAD, ZETA_SENS[el] * HA2CM, h.zeta / ZETA_SENS[el], h.F2 * HA2CM)
        FF[el].zeta === nothing || @printf(io, "  - zeta de Froese Fischer: %s cm^-1\n", FF[el].zeta)
    end
    println(io, "\n## zeta que pide el NIST con una sola zeta (Breit-Pauli p^2 sin CI)\n")
    for el in ELEMENTOS
        a, b = zeta_nist(el; excluir = 2), zeta_nist(el; excluir = 3)
        n = NIST_NIVELES[el]
        @printf(io, "- %s: sin 3P_1 -> zeta = %.1f (predice 3P_1 = %.1f, NIST %.1f); sin 3P_2 -> zeta = %.1f (predice 3P_2 = %.1f, NIST %.1f)\n",
                INFO[el].etiqueta, a.zeta, a.prediccion, n[2], b.zeta, b.prediccion, n[3])
    end
    println(io, "\n## zeta con termino tensorial dentro del 3P (ajuste exacto de los cuatro niveles)\n")
    for el in ELEMENTOS
        t = zeta_nist_tensorial(el)
        @printf(io, "- %s: zeta = %.1f, D = %.2f\n", INFO[el].etiqueta, t.zeta, t.D)
    end
    println(io, "\n## Escalamiento con Z: ajuste de ley de potencias y = A Z^p\n")
    println(io, "Minimos cuadrados de log y frente a log Z con los cuatro elementos (Z = 6, 14, 32, 50), ",
            "y exponente local entre elementos consecutivos. zeta y F^2 en cm^-1, <r^-3> en a0^-3.\n")
    let Zs = [INFO[el].Z for el in ELEMENTOS]
        series = ["zeta calculado"  => [HF[el].zeta * HA2CM for el in ELEMENTOS],
                  "zeta desnudo"    => [ALFA^2 / 2 * INFO[el].Z * HF[el].r3 * HA2CM for el in ELEMENTOS],
                  "zeta que pide el NIST" => [zeta_nist(el).zeta for el in ELEMENTOS],
                  "<r^-3>"          => [HF[el].r3 for el in ELEMENTOS],
                  "F^2(np,np)"      => [HF[el].F2 * HA2CM for el in ELEMENTOS]]
        for (nombre, ys) in series
            a = ajuste_potencia_Z(Zs, ys)
            loc = exponentes_locales(Zs, ys)
            @printf(io, "- %s: p = %.2f (R^2 = %.4f); local %s\n", nombre, a.p, a.R2,
                    join([@sprintf("%s->%s %.2f", INFO[ELEMENTOS[i]].etiqueta,
                                   INFO[ELEMENTOS[i + 1]].etiqueta, loc[i]) for i in eachindex(loc)], ", "))
        end
    end

    println(io, "\n## CI de valencia (alpha_d = 0, m de produccion)\n")
    for el in ELEMENTOS
        e = ci_de(el, ESTADO, 0.0)
        e === nothing && continue
        n = e["niveles_cm"]
        @printf(io, "- %s %s m = %d: E_corr = %.6e Ha (%.1f cm^-1), g(3P_2) con g_s = 2: %.6f, con g_s real: %.6f, niveles = %s\n",
                INFO[el].etiqueta, ESTADO, e["m"], e["E_corr_Ha"], e["E_corr_Ha"] * HA2CM, e["g_eff_3P2"], g_real(e),
                join([@sprintf("%.2f", x) for x in n], " / "))
    end
    println(io, "\n## Sensibilidad a la convencion de orbitales (alpha_d = 0, m = $(m_sens()))\n")
    for el in ELEMENTOS
        e3 = ci_de(el, ESTADO, 0.0; m = m_sens())
        ea = ci_de(el, ESTADO_SENSIBILIDAD, 0.0)
        (e3 === nothing || ea === nothing) && continue
        dif = [@sprintf("%s %+.2f%%", TERM_TXT[t], 100 * (ea["niveles_cm"][t] - e3["niveles_cm"][t]) / e3["niveles_cm"][t])
               for t in 2:5]
        @printf(io, "- %s: %s frente a %s: %s\n", INFO[el].etiqueta, ESTADO_SENSIBILIDAD, ESTADO, join(dif, ", "))
    end
    println(io, "\n## Estabilidad de los 3P_J con el espacio activo (3P, alpha_d = 0)\n")
    println(io, "Recorrido total de cada intervalo entre el primer y el ultimo m de la curva, y ",
            "excursion maxima respecto al valor en el m de produccion. En cm^-1.\n")
    for el in ELEMENTOS
        c = curva_convergencia(el)
        length(c) >= 2 || continue
        ms = [e["m"] for e in c]
        for (J, i) in ((1, 2), (2, 3))
            v = [e["niveles_cm"][i] for e in c]
            prod = v[findfirst(==(m_tabla(el)), ms)]
            @printf(io, "- %s 3P_%d: m = %d -> %d, %.2f -> %.2f (recorrido %.2f); excursion max. desde el m de produccion %.2f\n",
                    INFO[el].etiqueta, J, ms[1], ms[end], v[1], v[end], v[end] - v[1],
                    maximum(abs.(v .- prod)))
        end
    end

    println(io, "\n## Convergencia de los singletes con el espacio activo (3P, alpha_d = 0)\n")
    println(io, "Con los tres ultimos tamanos de cada curva: extrapolacion geometrica (incrementos en razon constante) ",
            "y de potencia (E(m) = E_inf + A m^-p, ajuste exacto). Una razon cercana a 1 hace inservible la geometrica.\n")
    for el in ELEMENTOS
        es = curva_convergencia(el)
        length(es) >= 3 || continue
        for (t, q) in ((4, "1D_2"), (5, "1S_0"))
            ms = [e["m"] for e in es[end-2:end]]
            ys = [e["niveles_cm"][t] for e in es[end-2:end]]
            g = extrap_geometrica(ys)
            p = extrap_potencia(ms, ys)
            nist = NIST_NIVELES[el][t]
            @printf(io, "- %s %s: m = %s -> %s cm^-1 | geometrica %s | potencia %s | NIST %.1f\n", INFO[el].etiqueta, q,
                    join(ms, ", "), join([@sprintf("%.1f", y) for y in ys], ", "),
                    g === nothing ? "no aplica" : @sprintf("%.1f (%+.1f%%, razon %.2f)", g.lim, 100 * (g.lim - nist) / nist, g.razon),
                    p === nothing ? "no aplica" : @sprintf("%.1f (%+.1f%%, p = %.1f)", p.lim, 100 * (p.lim - nist) / nist, p.p),
                    nist)
        end
    end
    println(io, "\n## Peso perturbativo de 1D_2 en el nivel 3P_2: w = zeta^2 / (2 Delta^2)\n")
    println(io, "Primer orden en el bloque J = 2 de Breit-Pauli, cuyo elemento fuera de la diagonal es ",
            "-zeta/sqrt(2); Delta = E(1D_2) - E(3P_2) entre los dos niveles de J = 2. Modelo: zeta calculado ",
            "y niveles del CI; NIST: zeta_NIST y niveles del ASD. 'diagonalizado' y 'del factor g' son los ",
            "pesos de la tabla de Lande, para comparar.\n")
    w_pert(z, n) = z^2 / (2 * (n[4] - n[3])^2)
    for el in ELEMENTOS
        e0 = ci_de(el, ESTADO, 0.0)
        e0 === nothing && continue
        zn, nn, gn = zeta_nist(el).zeta, NIST_NIVELES[el], NIST_LANDE_3P2[el]
        @printf(io, "- %s: modelo w = %.2f%% (diagonalizado %.2f%%); NIST w = %.2f%% (del factor g %s)",
                INFO[el].etiqueta, 100 * w_pert(HF[el].zeta * HA2CM, e0["niveles_cm"]), 100 * peso_1D(g_real(e0)),
                100 * w_pert(zn, nn), gn === nothing ? "sin dato" : @sprintf("%.2f%%", 100 * peso_1D(gn)))
        if haskey(BARRIDO_ALPHA, el)
            sel = alpha_elegido(el)
            if sel !== nothing
                ev = sel[2]
                @printf(io, "; con V_pol (alpha_d = %.2f) w = %.2f%% (diagonalizado %.2f%%)", sel[1],
                        100 * w_pert(ev["zeta_Ha"] * HA2CM, ev["niveles_cm"]), 100 * peso_1D(g_real(ev)))
            end
        end
        println(io)
    end
    println(io, "\n## Barridos de V_pol: error frente al NIST\n")
    for el in ELEMENTOS
        haskey(BARRIDO_ALPHA, el) || continue
        sel = alpha_elegido(el)
        for a in BARRIDO_ALPHA[el]
            e = ci_de(el, ESTADO, a)
            e === nothing && continue
            er = 100 .* error_3P(e["niveles_cm"], el)
            @printf(io, "- %s alpha_d = %.2f: 3P_1 %+.1f%%, 3P_2 %+.1f%%, 1D_2 = %.1f, 1S_0 = %.1f%s\n", INFO[el].etiqueta, a,
                    er[1], er[2], e["niveles_cm"][4], e["niveles_cm"][5], (sel !== nothing && sel[1] == a) ? "  <- elegido" : "")
        end
    end
    println(io, "\n## Cruce de alpha_d: el valor que iguala cada intervalo 3P_J al NIST\n")
    println(io, "Interpolacion lineal del barrido entre los dos alpha_d que acotan el nivel medido. ",
            "Un solo alpha_d basta para los dos intervalos si ambos cruces coinciden.\n")
    for el in ELEMENTOS
        haskey(BARRIDO_ALPHA, el) || continue
        as = vcat(0.0, BARRIDO_ALPHA[el])
        es = [ci_de(el, ESTADO, a) for a in as]
        ok = [i for i in eachindex(as) if es[i] !== nothing]
        for (J, k) in ((1, 2), (2, 3))
            meta = NIST_NIVELES[el][k]
            cruce = nothing
            for j in 1:length(ok) - 1
                i1, i2 = ok[j], ok[j + 1]
                v1, v2 = es[i1]["niveles_cm"][k], es[i2]["niveles_cm"][k]
                if (v1 - meta) * (v2 - meta) <= 0 && v2 != v1
                    cruce = as[i1] + (as[i2] - as[i1]) * (meta - v1) / (v2 - v1)
                    break
                end
            end
            @printf(io, "- %s 3P_%d (NIST %.1f cm^-1): %s\n", INFO[el].etiqueta, J, meta,
                    cruce === nothing ? "sin cruce dentro del barrido" : @sprintf("alpha_d = %.2f", cruce))
        end
    end

    println(io, "\n## Casi-degeneracion ns^2 np^2 <-> np^4 (dos configuraciones, orbitales congelados)\n")
    println(io, "G1(ns,np) y Delta en Ha; el resto en cm^-1. 'baja S-P' es cuanto desciende la separacion ",
            "1S - 3P por la mezcla con np^4. 'residuo' es el 1S_0 del CI de pareja menos el NIST, y ",
            "'sumado' lo que quedaria si la baja se sumara al CI: sobrecorrige, por eso no son aditivos.\n")
    for el in ELEMENTOS
        for a in vcat(0.0, haskey(BARRIDO_ALPHA, el) ? [alpha_elegido(el)[1]] : Float64[])
            d = casi_degeneracion(el, ESTADO, a)
            e = ci_de(el, ESTADO, a)
            e === nothing && continue
            res = e["niveles_cm"][5] - NIST_NIVELES[el][5]
            @printf(io, "- %s%s: G1 = %.5f, Delta = %.4f, peso de np^4 %.2f%% en 3P y %.2f%% en 1S, baja S-P %.0f, residuo 1S_0 %+.0f -> sumado %+.0f\n",
                    INFO[el].etiqueta, a == 0.0 ? "" : @sprintf(" (alpha_d = %.2f)", a),
                    d.G1, d.Delta, 100 * d.peso_P, 100 * d.peso_S, d.baja_S * HA2CM,
                    res, res - d.baja_S * HA2CM)
        end
    end

    println(io, "\n## Espacio de orbitales del CI: estado espurio del canal s y solapamiento con el core\n")
    println(io, "build_orbital_pool diagonaliza el Fock de core congelado sobre todos los splines, incluido el ",
            "primero, que no se anula en r = 0 (el SCF lo excluye). 'espurio' es el autovalor mas bajo del canal ",
            "s que descarta el filtro -Z^2; 'sin 1,N' es el mas bajo sin el primer y el ultimo spline. Despues, ",
            "max |<v|c>| y el peso maximo del core en un virtual del espacio activo (m de produccion), por l. ",
            "Detalle en tesis/diagnostico_estados_espurios.jl.\n")
    for el in ELEMENTOS
        ws = workspace(el)
        orbs = load(ruta_np2(archivo_scf(el, ESTADO, 0.0)))["orbitals"]
        core, n, Z = orbs[1:end-1], ws.basis.num_splines, INFO[el].Z
        ev1, ev2 = espectro_canal(ws, core, 0, 1:n), espectro_canal(ws, core, 0, 2:n-1)
        sol = solapamiento_con_core(ws, orbs, Z, m_tabla(el))
        @printf(io, "- %s: espurio %.1f Ha (filtro -Z^2 = %.0f, descartados %d); sin 1,N %.2f Ha; %s\n", INFO[el].etiqueta,
                minimum(ev1), -Z^2, count(<(-Z^2), ev1), minimum(ev2),
                join([@sprintf("l=%d: max |<v|c>| %.1e, peso %.1e", l, sol[l]...) for l in sort(collect(keys(sol)))], "; "))
    end

    println(io, "\n## C-DIIS: residuales de las trazas registradas (estado 3P, |dE| < 1e-10 Ha)\n")
    println(io, "De diis_trace_<elemento>_con_diis.csv. 'crudo' es max |FDS - SDF| y 'proyectado' el mismo ",
            "conmutador fuera del espacio ocupado, ambos en la ultima iteracion. 'piso desde' es la primera ",
            "iteracion en que el proyectado queda a menos del 10 % de su valor final; si coincide con la ",
            "ultima, el residual seguia bajando y no hay piso.\n")
    for el in ELEMENTOS
        r = residuales_diis(el)
        r === nothing && (push!(FALTAS, "traza con C-DIIS de $(INFO[el].prefijo)"); continue)
        @printf(io, "- %s: %d iteraciones; crudo %s, proyectado %s (crudo/proyectado = %.0f); piso desde la iteracion %d\n",
                INFO[el].etiqueta, r.iters, cient(r.crudo), cient(r.proy), r.crudo / r.proy, r.iter_piso)
    end

    println(io, "\n## Compresion del espectro por el CI de pareja (alpha_d = 0, m de produccion)\n")
    println(io, "Separaciones entre terminos del CI, sin espin-orbita, frente a las de Hartree-Fock ",
            "(6/25) F^2 y (15/25) F^2. En cm^-1.\n")
    for el in ELEMENTOS
        e = ci_de(el, ESTADO, 0.0)
        e === nothing && continue
        dD, dS = (e["E_1D_Ha"] - e["E_3P_Ha"]) * HA2CM, (e["E_1S_Ha"] - e["E_3P_Ha"]) * HA2CM
        F2 = e["F2_Ha"] * HA2CM
        @printf(io, "- %s: 1D - 3P = %.0f (HF %.0f, razon %.2f); 1S - 3P = %.0f (HF %.0f, razon %.2f, baja %.0f)\n",
                INFO[el].etiqueta, dD, 6F2 / 25, dD / (6F2 / 25), dS, 15F2 / 25, dS / (15F2 / 25), 15F2 / 25 - dS)
    end

    return String(take!(io))
end

function main_etapa3()
    mkpath(SALIDA)
    HF = Dict(el => propiedades_hf(el) for el in ELEMENTOS)
    ZETA_SENS = Dict(el => begin
                         d = load(ruta_np2(archivo_scf(el, ESTADO_SENSIBILIDAD)))
                         compute_zeta(d["R_grid"], d["V_eff"], d["P_$(INFO[el].n_val)p"])
                     end for el in ELEMENTOS)

    escribir_tabla("energia_global.tex", tabla_energia(HF)...; dir = SALIDA)
    escribir_tabla("parametros_radiales.tex", tabla_radiales(HF)...; dir = SALIDA)
    escribir_tabla("ionizacion.tex", tabla_ionizacion(HF)...; dir = SALIDA)
    escribir_tabla("niveles.tex", tabla_niveles(ELEMENTOS, false)...; dir = SALIDA)
    escribir_tabla("niveles_vpol.tex", tabla_niveles([el for el in ELEMENTOS if haskey(BARRIDO_ALPHA, el)], true)...;
                   dir = SALIDA)
    escribir_tabla("barrido_vpol.tex", tabla_barrido()...; dir = SALIDA)
    conv, costo, proc = tabla_convergencia()
    escribir_tabla("convergencia_ci.tex", conv, proc; dir = SALIDA)
    escribir_tabla("convergencia_ci_costo.tex", costo, proc; dir = SALIDA)
    escribir_tabla("convergencia_singletes.tex", tabla_convergencia_singletes()...; dir = SALIDA)
    escribir_tabla("cdiis.tex", tabla_cdiis()...; dir = SALIDA)
    escribir_tabla("parametros_base.tex", tabla_parametros_base()...; dir = SALIDA)
    escribir_tabla("zeta.tex", tabla_zeta(HF, ZETA_SENS)...; dir = SALIDA)
    escribir_tabla("acoplamiento.tex", tabla_acoplamiento(HF)...; dir = SALIDA)
    escribir_tabla("lande.tex", tabla_lande()...; dir = SALIDA)
    write(joinpath(SALIDA, "valores_texto.md"), valores_texto(HF, ZETA_SENS))
    println("  valores -> ", joinpath(SALIDA, "valores_texto.md"))

    if isempty(FALTAS)
        println("\nEtapa 3 completa: ninguna celda quedo sin dato.")
    else
        println("\nCeldas en \"--\" por falta de:")
        foreach(f -> println("  - ", f), sort(collect(FALTAS)))
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_etapa3()
end
