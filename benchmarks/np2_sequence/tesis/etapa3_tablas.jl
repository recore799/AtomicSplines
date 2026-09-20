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
    println(io, raw"\begin{tabular}{l l S[table-format=-4.8] S[table-format=-4.8] S[table-format=1.1e-1]}")
    println(io, raw"    \toprule")
    linea(io, raw"\textbf{Sistema}", raw"\textbf{Propiedad}", raw"{\textbf{Valor calculado}}",
          raw"{\textbf{Referencia F.F.}}", raw"{$\Delta$}")
    for el in ELEMENTOS
        h, ref = HF[el], FF[el]
        println(io, raw"    \midrule")
        println(io, "    \\multirow{3}{*}{$(etiqueta_Z(el))}")
        for (nombre, x, r) in ((raw"Energía total $E$", h.E, ref.E), (raw"Cinética $T$", h.T, ref.T),
                               ("Cociente virial", h.virial, ref.virial))
            d = r === nothing ? "{--}" : cient(abs(x - parse(Float64, r)))
            linea(io, "", nombre, @sprintf("%.8f", x), something(r, "{--}"), d)
        end
        push!(proc, proc_scf(h.archivo) * (ENERGIA_HF == :rayleigh ? ", E por cociente de Rayleigh" : ""))
    end
    println(io, raw"    \bottomrule", "\n", raw"\end{tabular}")
    push!(proc, "Referencias: Froese Fischer (1977), copiadas en tesis/config.jl.")
    return String(take!(io)), proc
end

function tabla_momentos(HF)
    io = IOBuffer(); proc = String[]
    println(io, raw"\begin{tabular}{l c S[table-format=1.6] S[table-format=1.6] S[table-format=1.1e-1]}")
    println(io, raw"    \toprule")
    linea(io, raw"\textbf{Átomo}", raw"\textbf{Capa}", raw"{\textbf{V.C. $\langle r^{-3} \rangle$}}",
          raw"{\textbf{F.F. $\langle r^{-3} \rangle$}}", raw"{$\Delta$}")
    println(io, raw"    \midrule")
    for el in ELEMENTOS
        h, r = HF[el], FF[el].r3
        d = r === nothing ? "{--}" : cient(abs(h.r3 - parse(Float64, r)))
        linea(io, etiqueta_Z(el), "\$$(INFO[el].n_val)p\$", @sprintf("%.6f", h.r3), something(r, "{--}"), d)
        push!(proc, proc_scf(h.archivo))
    end
    println(io, raw"    \bottomrule", "\n", raw"\end{tabular}")
    return String(take!(io)), proc
end

function tabla_slater(HF)
    io = IOBuffer(); proc = String[]
    println(io, raw"\begin{tabular}{l l S[table-format=1.8] S[table-format=1.8] S[table-format=1.1e-1]}")
    println(io, raw"    \toprule")
    linea(io, raw"\textbf{Átomo}", raw"\textbf{Multipolo}", raw"{\textbf{V.C. ($E_h$)}}",
          raw"{\textbf{F.F. ($E_h$)}}", raw"{$\Delta$}")
    for el in ELEMENTOS
        h, n = HF[el], INFO[el].n_val
        println(io, raw"    \midrule")
        println(io, "    \\multirow{2}{*}{$(INFO[el].etiqueta)}")
        for (k, x, r) in ((0, h.F0, FF[el].F0), (2, h.F2, FF[el].F2))
            d = r === nothing ? "{--}" : cient(abs(x - parse(Float64, r)))
            linea(io, "", "\$F^$k($(n)p, $(n)p)\$", @sprintf("%.8f", x), something(r, "{--}"), d)
        end
        push!(proc, proc_scf(h.archivo))
    end
    println(io, raw"    \bottomrule", "\n", raw"\end{tabular}")
    return String(take!(io)), proc
end

function tabla_koopmans(HF)
    io = IOBuffer(); proc = String[]
    println(io, raw"\begin{tabular}{lcccc}")
    println(io, raw"    \toprule")
    linea(io, raw"\textbf{Átomo}", raw"\textbf{SCF ($-\epsilon_{np}$, eV)}", raw"\textbf{CI ($\Delta E$, eV)}",
          raw"\textbf{NIST (eV)}", raw"\textbf{Error CI}")
    println(io, raw"    \midrule")
    for el in ELEMENTOS
        h, e = HF[el], ci_de(el, ESTADO, 0.0)
        scf = -h.eps * HA2EV
        ci = e === nothing ? nothing : (-h.eps - e["E_corr_Ha"]) * HA2EV
        nist = NIST_IONIZACION[el]
        err = ci === nothing ? "--" : @sprintf("%.1f\\%%", 100 * abs(ci - nist) / nist)
        linea(io, "$(INFO[el].etiqueta) (\$$(INFO[el].n_val)p\$)", fmt(scf, 2), fmt(ci, 2), fmt(nist, 2), err)
        push!(proc, proc_scf(h.archivo))
        e === nothing || push!(proc, proc_ci(e))
    end
    println(io, raw"    \bottomrule", "\n", raw"\end{tabular}")
    push!(proc, "CI (Delta E) = -eps_np - E_corr: el ion np^1 no tiene correlacion de pareja de valencia.")
    return String(take!(io)), proc
end

function tabla_niveles(els, con_vpol::Bool)
    io = IOBuffer(); proc = String[]
    cols = con_vpol ? 3 : 2
    println(io, "\\begin{tabular}{l|", join(fill("r"^cols, length(els)), "|"), "}")
    println(io, raw"    \toprule")
    enc = ["\\multicolumn{$cols}{c$(i < length(els) ? "|" : "")}{\\textbf{$(INFO[el].etiqueta) ($(el) I)}}"
           for (i, el) in enumerate(els)]
    linea(io, "", enc...)
    sub = con_vpol ? [raw"\textbf{HF+CI}", raw"\textbf{CI+$V_{\text{pol}}$}", raw"\textbf{NIST}"] :
                     [raw"\textbf{HF+CI}", raw"\textbf{NIST}"]
    linea(io, raw"\textbf{Término}", repeat(sub, length(els))...)
    println(io, raw"    \midrule")
    datos = Dict{String,Any}()
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
        datos[el] = (e0, ev)
    end
    for (t, nombre) in enumerate(TERMINOS)
        celdas = String[]
        for el in els
            e0, ev = datos[el]
            push!(celdas, e0 === nothing ? "--" : fmt(e0["niveles_cm"][t], 1))
            con_vpol && push!(celdas, ev === nothing ? "--" : fmt(ev["niveles_cm"][t], 1))
            push!(celdas, fmt(NIST_NIVELES[el][t], 1))
        end
        linea(io, nombre, celdas...)
    end
    println(io, raw"    \bottomrule", "\n", raw"\end{tabular}")
    return String(take!(io)), proc
end

function tabla_barrido(el)
    io = IOBuffer(); proc = String[]
    sel = alpha_elegido(el)
    println(io, raw"\begin{tabular}{ccccc}")
    println(io, raw"    \toprule")
    linea(io, raw"$\alpha_d$", raw"$\zeta_{np}$ (Ha)", raw"$^3P_1$ (cm$^{-1}$)", raw"$^3P_2$ (cm$^{-1}$)",
          CRITERIO_ALPHA == :solo_3P2 ? raw"Error $^3P_2$" : raw"Error rms $^3P_J$")
    println(io, raw"    \midrule")
    for a in vcat(0.0, get(BARRIDO_ALPHA, el, Float64[]))
        e = ci_de(el, ESTADO, a)
        marca = (sel !== nothing && sel[1] == a) ? raw"$^{*}$" : ""
        if e === nothing
            linea(io, fmt(a, 2) * marca, "--", "--", "--", "--")
        else
            n = e["niveles_cm"]
            linea(io, fmt(a, 2) * marca, fmt(e["zeta_Ha"], 5), fmt(n[2], 1), fmt(n[3], 1),
                  @sprintf("%.1f\\%%", 100 * costo_alpha(n, el)))
            push!(proc, proc_ci(e))
        end
    end
    println(io, raw"    \bottomrule", "\n", raw"\end{tabular}")
    push!(proc, "El asterisco marca el alpha_d elegido por el criterio $(CRITERIO_ALPHA); r_c = $(R_C) a0.")
    return String(take!(io)), proc
end

function tabla_ionizacion(el, HF)
    io = IOBuffer(); proc = String[]
    nist = NIST_IONIZACION[el]
    println(io, raw"\begin{tabular}{lccc}")
    println(io, raw"    \toprule")
    linea(io, "\\textbf{Modelo para $(el) ($(INFO[el].n_val)p)}", raw"\textbf{Ionización (eV)}",
          raw"\textbf{NIST (eV)}", raw"\textbf{Error}")
    println(io, raw"    \midrule")
    err(x) = x === nothing ? "--" : @sprintf("%.1f\\%%", 100 * abs(x - nist) / nist)
    h = HF[el]
    e0 = ci_de(el, ESTADO, 0.0)
    ip_hf = -h.eps * HA2EV
    ip_ci = e0 === nothing ? nothing : (-h.eps - e0["E_corr_Ha"]) * HA2EV
    linea(io, "HF (Koopmans)", fmt(ip_hf, 2), fmt(nist, 2), err(ip_hf))
    linea(io, "HF + CI (valencia)", fmt(ip_ci, 2), fmt(nist, 2), err(ip_ci))
    push!(proc, proc_scf(h.archivo))
    e0 === nothing || push!(proc, proc_ci(e0))
    sel = alpha_elegido(el)
    if sel !== nothing
        a, ev = sel
        eps_a = load(ruta_np2(ev["archivo"]))["orbitals"][end].energy
        ip_vp = (-eps_a - ev["E_corr_Ha"]) * HA2EV
        linea(io, @sprintf("HF + CI + \$V_{\\text{pol}}\$ (\$\\alpha_d = %.2f\$)", a), fmt(ip_vp, 2), fmt(nist, 2), err(ip_vp))
        push!(proc, proc_scf(ev["archivo"]), proc_ci(ev))
    else
        linea(io, raw"HF + CI + $V_{\text{pol}}$", "--", fmt(nist, 2), "--")
    end
    println(io, raw"    \bottomrule", "\n", raw"\end{tabular}")
    return String(take!(io)), proc
end

function tabla_convergencia()
    io = IOBuffer(); io2 = IOBuffer(); proc = String[]
    por_el = Dict(el => Dict(e["m"] => e for e in CI if e["archivo"] == archivo_scf(el, ESTADO, 0.0) &&
                                                      e["sha256_entrada"] == sha256_de(ruta_np2(archivo_scf(el, ESTADO, 0.0))))
                  for el in ELEMENTOS)
    ms = sort(unique(vcat([collect(keys(d)) for d in values(por_el)]...)))
    println(io, raw"\begin{tabular}{rrr", "r"^length(ELEMENTOS), "}")
    println(io, raw"    \toprule")
    linea(io, raw"$m$", raw"Orbitales", raw"CSF $^3P$", ["\$E_{\\text{corr}}\$ $(el) (mHa)" for el in ELEMENTOS]...)
    println(io, raw"    \midrule")
    println(io2, raw"\begin{tabular}{r", "r"^(2 * length(ELEMENTOS)), "}")
    println(io2, raw"    \toprule")
    linea(io2, raw"$m$", ["\$R^k\$ $(el)" for el in ELEMENTOS]..., ["\$t\$ $(el) (s)" for el in ELEMENTOS]...)
    println(io2, raw"    \midrule")
    for m in ms
        ref = first(d[m] for d in values(por_el) if haskey(d, m))
        linea(io, string(m), string(ref["n_orb"]), string(ref["n_csf"][1]),
              [haskey(por_el[el], m) ? fmt(1e3 * por_el[el][m]["E_corr_Ha"], 3) : "--" for el in ELEMENTOS]...)
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
    linea(io, raw"\textbf{Elemento}", raw"\textbf{Iter. sin C-DIIS}", raw"\textbf{Iter. con C-DIIS}",
          raw"\textbf{Factor}", raw"$|\Delta E|$ entre ambos (Ha)")
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
              dE < 10.0^(-d) ? "\$<10^{-$d}\$" : cient(dE))
        push!(proc, "$(s.archivo), $(c.archivo)")
    end
    println(io, raw"    \bottomrule", "\n", raw"\end{tabular}")
    push!(proc, "Trazas de diis_benchmark.jl (estado 3P, |dE| < 1e-10 Ha, mismo punto de partida). " *
                "Las sin C-DIIS de Ge y Sn se reconstruyeron del log impreso, con E a 8 decimales.")
    return String(take!(io)), proc
end

function tabla_zeta(HF, ZETA_SENS)
    io = IOBuffer(); proc = String[]
    println(io, raw"\begin{tabular}{lrccccccc}")
    println(io, raw"    \toprule")
    linea(io, raw"\textbf{Elemento}", raw"$Z$", raw"$(\alpha Z)^2$", raw"$\zeta_{\text{desnudo}}$",
          raw"$\zeta_{\text{calc}}$", raw"$\zeta_{\text{NIST}}$", raw"$\zeta_{\text{calc}}/\zeta_{\text{NIST}}$",
          raw"$R_{\text{NIST}}$", raw"$R_{\text{modelo}}$")
    println(io, raw"    \midrule")
    for el in ELEMENTOS
        h, Z = HF[el], INFO[el].Z
        zb = ALFA^2 / 2 * Z * h.r3 * HA2CM
        zc = h.zeta * HA2CM
        zn = zeta_nist(el).zeta
        n = NIST_NIVELES[el]
        e0 = ci_de(el, ESTADO, 0.0)
        Rmod = e0 === nothing ? nothing : (e0["niveles_cm"][3] - e0["niveles_cm"][2]) / e0["niveles_cm"][2]
        linea(io, INFO[el].etiqueta, string(Int(Z)), fmt((ALFA * Z)^2, 4), fmt(zb, 1), fmt(zc, 1), fmt(zn, 1),
              fmt(zc / zn, 3), fmt((n[3] - n[2]) / n[2], 3), fmt(Rmod, 3))
        push!(proc, proc_scf(h.archivo))
        e0 === nothing || push!(proc, proc_ci(e0))
    end
    println(io, raw"    \bottomrule", "\n", raw"\end{tabular}")
    push!(proc, "zeta en cm^-1. zeta_desnudo = (alpha^2/2) Z <r^-3>. zeta_NIST: la zeta unica que reproduce " *
                "3P_2, 1D_2 y 1S_0 del NIST en Breit-Pauli de p^2 sin CI. R = (3P_2 - 3P_1)/(3P_1 - 3P_0); " *
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
    pct(x) = x === nothing ? "--" : @sprintf("%.2f\\%%", 100 * x)
    println(io, raw"\begin{tabular}{lcccccc}")
    println(io, raw"    \toprule")
    linea(io, raw"\textbf{Elemento}", raw"$g$ HF+CI", raw"$g$ CI+$V_{\text{pol}}$", raw"$g$ NIST",
          raw"$^1D_2$ HF+CI", raw"$^1D_2$ CI+$V_{\text{pol}}$", raw"$^1D_2$ NIST")
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
    escribir_tabla("momentos_inversos.tex", tabla_momentos(HF)...; dir = SALIDA)
    escribir_tabla("integrales_slater.tex", tabla_slater(HF)...; dir = SALIDA)
    escribir_tabla("koopmans.tex", tabla_koopmans(HF)...; dir = SALIDA)
    escribir_tabla("niveles_ligeros.tex", tabla_niveles(["C", "Si"], false)...; dir = SALIDA)
    escribir_tabla("niveles_pesados.tex", tabla_niveles(["Ge", "Sn"], true)...; dir = SALIDA)
    for el in ELEMENTOS
        haskey(BARRIDO_ALPHA, el) || continue
        escribir_tabla("barrido_vpol_$(el).tex", tabla_barrido(el)...; dir = SALIDA)
        escribir_tabla("ionizacion_$(el).tex", tabla_ionizacion(el, HF)...; dir = SALIDA)
    end
    conv, costo, proc = tabla_convergencia()
    escribir_tabla("convergencia_ci.tex", conv, proc; dir = SALIDA)
    escribir_tabla("convergencia_ci_costo.tex", costo, proc; dir = SALIDA)
    escribir_tabla("convergencia_singletes.tex", tabla_convergencia_singletes()...; dir = SALIDA)
    escribir_tabla("cdiis.tex", tabla_cdiis()...; dir = SALIDA)
    escribir_tabla("zeta.tex", tabla_zeta(HF, ZETA_SENS)...; dir = SALIDA)
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
