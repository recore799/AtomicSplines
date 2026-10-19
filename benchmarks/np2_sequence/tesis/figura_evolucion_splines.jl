# =============================================================================
#  Figura 4.1 del manuscrito: B-splines de grado 0 a 3 por la recurrencia de Cox-de Boor
#
#  Ilustracion, no resultado: no lee ningun .jld2 ni toca el manifiesto. Evalua B_i^p(r) con la
#  recurrencia de la Seccion 4.1.1, tal como la escribe el texto, sobre nudos uniformes t_i = i,
#  y dibuja un panel por grado. La version anterior de la figura, hecha a mano para un cartel,
#  no tenia generador; esta conserva su contenido y usa la notacion del texto (r y nudos t_i).
#
#  Uso: julia --project=. benchmarks/np2_sequence/tesis/figura_evolucion_splines.jl
#       ... --salida /ruta/dir      (para pruebas)
# =============================================================================

ENV["GKSwstype"] = "100"            # GR sin abrir ventanas
using Plots

const REPO_FIG = normpath(joinpath(@__DIR__, "..", "..", ".."))

function argumento_fig(nombre, defecto)
    i = findfirst(==(nombre), ARGS)
    return (i === nothing || i == length(ARGS)) ? defecto : ARGS[i + 1]
end

const SALIDA_FIG = argumento_fig("--salida", joinpath(REPO_FIG, "docs", "figures"))
const NUDOS      = collect(0.0:9.0)            # t_i = i
const VENTANA    = (1.5, 6.5)                  # tramo de la malla que se muestra
const TITULOS    = ["a) Grado p = 0 (escalón)", "b) Grado p = 1 (lineal)",
                    "c) Grado p = 2 (cuadrático)", "d) Grado p = 3 (cúbico)"]

"""B_i^p(r) por la recurrencia de Cox-de Boor; `i` es el indice del nudo t_i (base 0)."""
function cox_de_boor(i, p, r, t = NUDOS)
    if p == 0
        return t[i + 1] <= r < t[i + 2] ? 1.0 : 0.0
    end
    izq = (r - t[i + 1]) / (t[i + p + 1] - t[i + 1]) * cox_de_boor(i, p - 1, r, t)
    der = (t[i + p + 2] - r) / (t[i + p + 2] - t[i + 2]) * cox_de_boor(i + 1, p - 1, r, t)
    return izq + der
end

function panel(p)
    r = collect(range(VENTANA[1], VENTANA[2], length = 1001))
    marcas = filter(t -> VENTANA[1] <= t <= VENTANA[2], NUDOS)
    pl = plot(title = TITULOS[p + 1], ylabel = "\$B_i^{$p}(r)\$", xlims = VENTANA, ylims = (0, 1.1),
              xticks = (marcas, ["\$t_{$(Int(t))}\$" for t in marcas]), legend = false,
              titlefontsize = 12, tickfontsize = 10, guidefontsize = 11, grid = false)
    for i in 0:(length(NUDOS) - p - 2)       # todos los B_i^p que caben en el vector de nudos
        y = [cox_de_boor(i, p, x) for x in r]
        any(>(0), y) || continue
        plot!(pl, r, y; lw = 1.5, fillrange = 0, fillalpha = 0.12, seriescolor = i + 1)
    end
    return pl
end

function main_figura_evolucion()
    paneles = [panel(p) for p in 0:3]
    xlabel!(paneles[end], "\$r\$ (las marcas son los nudos \$t_i\$)")
    mkpath(SALIDA_FIG)
    ruta = joinpath(SALIDA_FIG, "figura_evolucion_splines.pdf")
    savefig(plot(paneles...; layout = (4, 1), size = (700, 950), margin = 3Plots.mm), ruta)
    println("  figura -> ", ruta)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_figura_evolucion()
end
