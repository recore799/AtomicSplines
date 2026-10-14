# =============================================================================
#  Figura 4.3 del manuscrito: distribucion uniforme de nudos frente a la malla exponencial
#
#  Ilustracion, no resultado: no lee ningun .jld2 ni toca el manifiesto. Dibuja la base de
#  B-splines que construye el motor (misma secuencia de nudos que generate_basis y misma
#  evaluacion que evaluate_orbital) con pocos intervalos, para que las funciones se distingan.
#
#  GAMMA = 4.0 es el valor de la malla del estanio en la tabla de parametros de la base. La
#  version anterior de la figura, hecha a mano con gamma = 2.5, no tenia generador.
#
#  Uso: julia --project=. benchmarks/np2_sequence/tesis/figura_comparacion_mallas.jl
#       ... --gamma 4.0 --salida /ruta/dir      (para pruebas)
# =============================================================================

using AtomicSplines
ENV["GKSwstype"] = "100"            # GR sin abrir ventanas
using Plots
using Printf

const REPO_FIG = normpath(joinpath(@__DIR__, "..", "..", ".."))

function argumento_fig(nombre, defecto)
    i = findfirst(==(nombre), ARGS)
    return (i === nothing || i == length(ARGS)) ? defecto : ARGS[i + 1]
end

const R_MAX_FIG  = 20.0
const N_INTERV   = 15
const ORDEN      = 4                # splines cubicos: orden k = 4, grado 3
const GAMMA      = parse(Float64, argumento_fig("--gamma", "4.0"))
const SALIDA_FIG = argumento_fig("--salida", joinpath(REPO_FIG, "docs", "figures"))

"""Base con nudos interiores dados y multiplicidad `ORDEN` en los extremos, como generate_basis."""
function base_con_nudos(r_nudos)
    nudos = vcat(fill(r_nudos[1], ORDEN - 1), r_nudos, fill(r_nudos[end], ORDEN - 1))
    return BSplineBasis{ORDEN}(length(nudos) - ORDEN, nudos, Float64[], Float64[])
end

nudos_uniformes() = collect(range(0.0, R_MAX_FIG, length = N_INTERV + 1))
nudos_exponenciales(γ) = unique(generate_basis(R_MAX_FIG, N_INTERV, Val(ORDEN); γ = γ).knots)

"""Puntos de muestreo: los mismos por intervalo, para resolver los splines estrechos del origen."""
muestreo(r_nudos; por_intervalo = 60) =
    unique(vcat([collect(range(r_nudos[i], r_nudos[i + 1], length = por_intervalo)) for i in 1:(length(r_nudos) - 1)]...))

function panel(r_nudos, color, titulo)
    base = base_con_nudos(r_nudos)
    r = muestreo(r_nudos)
    p = plot(title = titulo, xlabel = "r (u.a.)", ylabel = "Bᵢ(r)", xlims = (0, R_MAX_FIG),
             ylims = (0, 1.05), legend = false, grid = false)
    for i in 1:base.num_splines
        c = zeros(base.num_splines)
        c[i] = 1.0
        plot!(p, r, evaluate_orbital(base, c, r); color = color, lw = 1.5, fillrange = 0, fillalpha = 0.08)
    end
    for t in r_nudos                # marcas de los nudos sobre el eje
        plot!(p, [t, t], [0.0, 0.05]; color = :black, lw = 1.5)
    end
    return p
end

function main_figura_mallas()
    a = panel(nudos_uniformes(), :blue, "a) Distribución uniforme")
    b = panel(nudos_exponenciales(GAMMA), :red, @sprintf("b) Distribución exponencial (γ = %g)", GAMMA))
    mkpath(SALIDA_FIG)
    ruta = joinpath(SALIDA_FIG, "figura_comparacion_mallas.pdf")
    savefig(plot(a, b; layout = (2, 1), size = (800, 600), margin = 4Plots.mm), ruta)
    println("  figura -> ", ruta)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_figura_mallas()
end
