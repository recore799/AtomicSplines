# =============================================================================
#  Etapa 4: figuras del capitulo de resultados (un par de minutos, casi todo es compilar
#  Plots la primera vez)
#
#  Regenera en docs/figures/ las siete figuras que usa resultados.tex, ahora desde los .jld2
#  registrados y resultados_ci.toml en lugar de constantes escritas a mano o de archivos
#  sueltos de examples/scratch/, y agrega cuatro nuevas: convergencia del CI, convergencia de
#  C-DIIS, razon de zeta y barrido de V_pol.
#
#  Uso: julia --project=. benchmarks/np2_sequence/tesis/etapa4_figuras.jl
#       ... --ci /ruta/prueba.toml --m 3 --salida /ruta/dir      (para pruebas)
# =============================================================================

include(joinpath(@__DIR__, "comun.jl"))
isdefined(Main, :run_full_ci) || include(joinpath(NP2, "np2_ci_full.jl"))
include(joinpath(@__DIR__, "analisis_zeta.jl"))
ENV["GKSwstype"] = "100"            # GR sin abrir ventanas
using Plots

const CI     = leer_ci(argumento("--ci", RESULTADOS_CI))
const M      = parse(Int, argumento("--m", string(M_PRODUCCION)))
const SALIDA = argumento("--salida", DIR_FIGURAS)
const MAN    = manifiesto()
const COLOR  = Dict("C" => :blue, "Si" => :green, "Ge" => :red, "Sn" => :orange)
const ESTILO = Dict("C" => :solid, "Si" => :dash, "Ge" => :dot, "Sn" => :dashdot)
const NOMBRES_NIVELES = ["³P₀", "³P₁", "³P₂", "¹D₂", "¹S₀"]

# Eje en sqrt(r) con las marcas de los scripts originales de examples/scratch/, tomadas de Froese
# Fischer: abre la region cercana al nucleo sin perder la cola de valencia.
const R_PLOT  = 12.5
const R_TICKS = [0.0, 0.125, 0.5, 1.125, 2.0, 3.125, 4.5, 6.125, 8.0, 10.125, 12.5]
const XTICKS  = (sqrt.(R_TICKS), ["0", "0.12", "0.5", "1.12", "2", "3.12", "4.5", "6.12", "8", "10.1", "12.5"])

etiqueta_Z(el) = "$(INFO[el].etiqueta) (Z=$(Int(INFO[el].Z)))"

function ci_fig(el, estado, a)
    archivo = archivo_scf(el, estado, a)
    return registrado(archivo, MAN) ? buscar_ci(CI, archivo, M) : nothing
end

function datos_hf(el)
    archivo = archivo_scf(el, ESTADO)
    registrado(archivo, MAN) || error("$archivo no esta registrado.")
    verificar_malla(el, archivo, MAN)
    return load(ruta_np2(archivo))
end

"""Signo de cada orbital: primer lobulo significativo positivo."""
function fase(P)
    i = findfirst(x -> abs(x) > 1e-4, P)
    return (i !== nothing && P[i] < 0) ? -P : P
end

"""
Guarda la figura. Los PDF no se versionan (`*.pdf` en .gitignore), asi que la primera vez que
se pisa una figura la version anterior se copia a `anteriores/`: es la unica copia que queda
y sirve para comparar.
"""
function guardar(p, nombre)
    mkpath(SALIDA)
    ruta = joinpath(SALIDA, nombre)
    previa = joinpath(SALIDA, "anteriores", nombre)
    if isfile(ruta) && !isfile(previa)
        mkpath(dirname(previa))
        cp(ruta, previa)
    end
    savefig(p, ruta)
    println("  figura -> ", ruta)
end

function figuras_orbitales(D)
    plog(V) = sign(V) * log10(1 + abs(V) / 0.1)
    yv = [0, -0.5, -1, -5, -10, -50, -100, -500, -1000, -3500]
    comun = (xlabel = "r (u.a.)", size = (850, 500), margin = 5Plots.mm, xticks = XTICKS,
             xlims = (0, sqrt(R_PLOT)))
    # El fondo del pozo de V_rad baja con Z. El limite inferior sale de los datos: el -300 fijo de
    # plot_np2_sequence.jl cortaba al estanio, cuyo minimo ronda -400 Ha.
    fondo = minimum(minimum(d["V_eff"][d["R_grid"] .<= R_PLOT] .+ 1.0 ./ d["R_grid"][d["R_grid"] .<= R_PLOT] .^ 2)
                    for d in values(D))
    p1 = plot(; comun..., title = "Funciones de Onda Radiales de Valencia (Secuencia np²)",
              ylabel = "P_np(r)", legend = :topright)
    p2 = plot(; comun..., title = "Potenciales Centrales Efectivos (Secuencia np²)",
              ylabel = "V_eff(r) (Ha)", legend = :bottomright,
              yticks = (plog.(yv), string.(yv)), ylims = (plog(-3600), plog(2)))
    p3 = plot(; comun..., title = "Potenciales Radiales Efectivos (V_eff + 1/r²)",
              ylabel = "V_rad(r) (Ha)", legend = :bottomright, ylims = (1.1 * fondo, 50))
    paneles_orb, paneles_dens = [], []
    for el in ELEMENTOS
        i = INFO[el]
        R = D[el]["R_grid"]
        idx = findall(<=(R_PLOT), R)
        x = sqrt.(R[idx])
        V = D[el]["V_eff"][idx]
        plot!(p1, x, fase(D[el]["P_$(i.n_val)p"][idx]); label = "$(i.etiqueta) $(i.n_val)p",
              lw = 2.5, color = COLOR[el], linestyle = ESTILO[el])
        plot!(p2, x, plog.(V); label = "$(i.etiqueta) V_eff (Z=$(Int(i.Z)))", lw = 2, color = COLOR[el])
        plot!(p3, x, V .+ 1.0 ./ R[idx] .^ 2; label = "$(i.etiqueta) V_rad (Z=$(Int(i.Z)))", lw = 2,
              color = COLOR[el])

        # Todos los orbitales con la base del propio elemento. La num_splines tiene que
        # coincidir: plot_all_orbitals.jl usaba N_elems = num_splines - K, un elemento de
        # menos, y evaluaba los coeficientes sobre una malla de nudos que no era la suya.
        basis = generate_basis(R_MAX, i.N, Val(i.K); γ = i.gamma)
        eje_x = el == last(ELEMENTOS) ? "r (u.a.)" : ""
        po = plot(title = etiqueta_Z(el), ylabel = "P_nl(r)", xlabel = eje_x, xticks = XTICKS,
                  xlims = (0, sqrt(R_PLOT)), legend = :topright, margin = 3Plots.mm, left_margin = 8Plots.mm)
        dens = zeros(length(idx))
        for o in D[el]["orbitals"]
            basis.num_splines == length(o.coeffs) ||
                error("$el: la base tiene $(basis.num_splines) splines y el orbital $(length(o.coeffs)).")
            P = evaluate_orbital(basis, o.coeffs, R[idx])
            plot!(po, x, fase(P); label = "$(o.n)$("spdf"[o.l + 1])", lw = 2)
            dens .+= o.occ .* P .^ 2
        end
        push!(paneles_orb, po)
        push!(paneles_dens, plot(x, dens; title = "Densidad Total $(etiqueta_Z(el))", ylabel = "D(r)", xlabel = eje_x,
                                 xticks = XTICKS, xlims = (0, sqrt(R_PLOT)), legend = false, lw = 2,
                                 color = :black, fill = (0, 0.2, COLOR[el]), margin = 3Plots.mm, left_margin = 8Plots.mm))
    end
    hline!(p3, [0]; color = :black, lw = 1, linestyle = :dash, label = "")
    guardar(p1, "np2_valence_wavefunctions.pdf")
    guardar(p2, "np2_effective_potentials.pdf")
    guardar(p3, "np2_radial_potentials.pdf")
    guardar(plot(paneles_orb...; layout = (4, 1), size = (850, 1300)), "np2_all_orbitals.pdf")
    guardar(plot(paneles_dens...; layout = (4, 1), size = (850, 1300)), "np2_total_densities.pdf")
end

function figuras_zeta(HFV)
    Zs = [INFO[el].Z for el in ELEMENTOS]
    xt = (Zs, string.(Int.(Zs)))
    zc = [HFV[el].zeta * HA2CM for el in ELEMENTOS]
    zn = [zeta_nist(el).zeta for el in ELEMENTOS]
    zb = [ALFA^2 / 2 * INFO[el].Z * HFV[el].r3 * HA2CM for el in ELEMENTOS]
    F2 = [HFV[el].F2 * HA2CM for el in ELEMENTOS]

    p = plot(Zs, zc; xscale = :log10, yscale = :log10, xticks = xt, marker = :circle, lw = 2,
             color = :red, label = "ζ_np calculado", xlabel = "Z", ylabel = "Energía (cm⁻¹)",
             legend = :bottomright, size = (850, 500), margin = 5Plots.mm)
    scatter!(p, Zs, zn; marker = :diamond, color = :black, label = "ζ que pide el NIST")
    plot!(p, Zs, F2; marker = :square, lw = 2, color = :blue, label = "F²(np,np)")
    guardar(p, "scaling_law.pdf")

    q = plot(Zs, zc ./ zn; xticks = xt, marker = :circle, lw = 2, color = :red,
             label = "ζ calculado / ζ NIST", xlabel = "Z", ylabel = "razón", legend = :topright,
             size = (850, 500), margin = 5Plots.mm)
    plot!(q, Zs, zb ./ zn; marker = :square, lw = 2, linestyle = :dash, color = :gray,
          label = "ζ de núcleo desnudo / ζ NIST")
    hline!(q, [1.0]; color = :black, linestyle = :dot, label = "")
    guardar(q, "zeta_razon.pdf")
end

function figura_niveles()
    paneles = []
    for el in ELEMENTOS
        e = ci_fig(el, ESTADO, 0.0)
        p = plot(title = etiqueta_Z(el), legend = false, xticks = ([0.5, 1.8], ["HF+CI", "NIST"]),
                 xlims = (-0.2, 2.9), ylabel = "Energía (cm⁻¹)")
        for (x0, niv, color) in ((0.0, e === nothing ? nothing : e["niveles_cm"], :blue),
                                 (1.3, NIST_NIVELES[el], :black))
            niv === nothing && continue
            for E in niv
                plot!(p, [x0, x0 + 1.0], [E, E]; lw = 2, color = color)
            end
        end
        for (t, E) in enumerate(NIST_NIVELES[el])
            t in (2, 3) && continue          # 3P_1 y 3P_2 se encimarian con 3P_0 a esta escala
            annotate!(p, 2.35, E, text(NOMBRES_NIVELES[t] * (t == 1 ? "," * NOMBRES_NIVELES[2][end:end] *
                                                               "," * NOMBRES_NIVELES[3][end:end] : ""), 8, :left))
        end
        push!(paneles, p)
    end
    guardar(plot(paneles...; layout = (1, 4), size = (1300, 500), margin = 5Plots.mm), "np2_energy_levels.pdf")
end

function figura_convergencia()
    p = plot(xlabel = "m (orbitales por cada l ≤ $(LMAX_CI))", ylabel = "E_corr (mHa)",
             legend = :topright, size = (850, 500), margin = 5Plots.mm)
    for el in ELEMENTOS
        archivo = archivo_scf(el, ESTADO)
        registrado(archivo, MAN) || continue
        sha = sha256_de(ruta_np2(archivo))
        es = sort([e for e in CI if e["archivo"] == archivo && e["sha256_entrada"] == sha]; by = e -> e["m"])
        isempty(es) && continue
        plot!(p, [e["m"] for e in es], [1e3 * e["E_corr_Ha"] for e in es]; marker = :circle, lw = 2,
              color = COLOR[el], label = INFO[el].etiqueta)
    end
    guardar(p, "convergencia_ci.pdf")
end

function figura_cdiis()
    paneles = []
    for el in ("Ge", "Sn")
        p = plot(title = INFO[el].etiqueta, xlabel = "iteración", yscale = :log10, legend = :topright)
        for (tag, color, nombre) in (("sin_diis", :gray, "sin C-DIIS"), ("con_diis", :red, "con C-DIIS"))
            ruta = joinpath(NP2, "diis_trace_$(INFO[el].prefijo)_$(tag).csv")
            isfile(ruta) || continue
            filas = [parse.(Float64, split(l, ",")) for l in readlines(ruta) if !isempty(strip(l)) && !startswith(l, "#") && !startswith(l, "iter")]
            it = [f[1] for f in filas]
            for (col, estilo, que) in ((3, :solid, "|ΔE| (Ha)"), (4, :dash, "residual proyectado"))
                y = [f[col] for f in filas]
                ok = isfinite.(y) .& (y .> 0)
                any(ok) && plot!(p, it[ok], y[ok]; lw = 2, color = color, linestyle = estilo,
                                 label = "$nombre, $que")
            end
        end
        push!(paneles, p)
    end
    guardar(plot(paneles...; layout = (1, 2), size = (1200, 450), margin = 5Plots.mm), "cdiis_convergencia.pdf")
end

function figura_barrido()
    paneles = []
    for el in ELEMENTOS
        haskey(BARRIDO_ALPHA, el) || continue
        as, e1, e2 = Float64[], Float64[], Float64[]
        for a in vcat(0.0, BARRIDO_ALPHA[el])
            e = ci_fig(el, ESTADO, a)
            e === nothing && continue
            n, ref = e["niveles_cm"], NIST_NIVELES[el]
            push!(as, a)
            push!(e1, 100 * (n[2] - ref[2]) / ref[2])
            push!(e2, 100 * (n[3] - ref[3]) / ref[3])
        end
        p = plot(title = INFO[el].etiqueta, xlabel = "α_d (u.a.)", ylabel = "error frente al NIST (%)",
                 legend = :topleft)
        if !isempty(as)
            plot!(p, as, e1; marker = :circle, lw = 2, label = "³P₁")
            plot!(p, as, e2; marker = :square, lw = 2, label = "³P₂")
        end
        hline!(p, [0.0]; color = :black, linestyle = :dot, label = "")
        push!(paneles, p)
    end
    guardar(plot(paneles...; layout = (1, length(paneles)), size = (1200, 450), margin = 5Plots.mm),
            "barrido_vpol.pdf")
end

function main_etapa4()
    D = Dict(el => datos_hf(el) for el in ELEMENTOS)
    HFV = Dict{String,Any}()
    for el in ELEMENTOS
        ws = workspace(el)
        v = D[el]["orbitals"][end]
        HFV[el] = (r3 = dot(v.coeffs, ws.R_inv3 * v.coeffs), F2 = compute_Rk(ws, v, v, v, v, 2),
                   zeta = compute_zeta(D[el]["R_grid"], D[el]["V_eff"], D[el]["P_$(INFO[el].n_val)p"]))
    end
    figuras_orbitales(D)
    figuras_zeta(HFV)
    figura_niveles()
    figura_convergencia()
    figura_cdiis()
    figura_barrido()
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_etapa4()
end
