# =============================================================================
#  Banco de pruebas de C-DIIS para la secuencia np^2
#
#  Corre cada elemento DOS veces con el mismo codigo, el mismo punto de partida y la
#  misma tolerancia, y lo unico que cambia entre las dos corridas es el acelerador:
#
#    sin DIIS : SCF simple, con level shift permanente de 3 Ha en Ge y Sn.
#    con DIIS : level shift solo hasta que max_canales||FDS-SDF||_inf < diis_thresh,
#               y de ahi en adelante extrapolacion de Pulay sin shift.
#
#  Corre con save = false: NO reescribe los *_rohf_results_*.jld2, que son entrada del
#  trabajo de CI. La energia de referencia se LEE de esos archivos, no se copia aqui.
#
#  Uso:  julia --project=../.. diis_benchmark.jl [Carbon] [Silicon] [Germanium] [Tin]
#        (sin argumentos corre los cuatro; Ge y Sn tardan minutos)
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using JLD2, Printf, AtomicSplines

include("carbon_rohf.jl")
include("silicon_rohf.jl")
include("germanium_rohf.jl")
include("tin_rohf.jl")

const TOL      = 1e-10
const R_MAX    = 30.0
const ESTADO   = "3P"

const CASES = Dict(
    "Carbon"    => (solve_carbon_rohf,    "carbon_rohf_results_3P_R30.0.jld2",     600),
    "Silicon"   => (solve_silicon_rohf,   "silicon_rohf_results_3P_R30.0.jld2",    600),
    "Germanium" => (solve_germanium_rohf, "germanium_rohf_results_3P_R30.0.jld2", 1000),
    "Tin"       => (solve_tin_rohf,       "tin_rohf_results_3P_R30.0.jld2",       1000),
)
const ORDER = ["Carbon", "Silicon", "Germanium", "Tin"]

function run_case(name::String; verbose::Bool = false)
    solver, reffile, max_iter = CASES[name]
    ref = joinpath(@__DIR__, reffile)
    E_ref = isfile(ref) ? load(ref, "E_total") : NaN

    rows = Dict{Bool,Any}()
    for use_diis in (false, true)
        @printf("\n##### %s -- %s #####\n", name, use_diis ? "con C-DIIS" : "sin C-DIIS")
        t0 = time()
        r = solver(R_MAX; estado = ESTADO, use_diis = use_diis, tol = TOL,
                   max_iter = max_iter, save = false, verbose = verbose)
        rows[use_diis] = (r = r, wall = time() - t0)
        @printf("  -> E = %.10f Ha | %d iteraciones | %.1f s%s\n",
                r.E_total, r.iters, time() - t0,
                r.converged ? "" : "  [NO CONVERGIO]")
    end
    return (name = name, E_ref = E_ref, off = rows[false], on = rows[true])
end

"""
    write_traces(results; dir)

Vuelca la traza de cada corrida a CSV para la figura de convergencia. Una fila por
iteracion con la energia, |dE| y los DOS residuales: el proyectado, que es el que
gobierna el cambio de regimen, y el conmutador crudo FDS-SDF, que se queda en un piso
no nulo. Esa diferencia es el argumento que hay que mostrar en la tesis.
"""
function write_traces(results; dir::String = ".")
    for res in results, (tag, row) in (("sin_diis", res.off), ("con_diis", res.on))
        path = joinpath(dir, "diis_trace_$(lowercase(res.name))_$(tag).csv")
        open(path, "w") do io
            println(io, "iter,E_total,delta_E,resid_proj,resid_raw")
            for (it, E, d, rp, rr) in row.r.trace
                @printf(io, "%d,%.12f,%.6e,%.6e,%.6e\n", Int(it), E, d, rp, rr)
            end
        end
        println("  traza -> ", path)
    end
end

function report(results)
    println("\n\n", "="^96)
    @printf(" C-DIIS en el ROHF de la secuencia np^2 -- estado %s, tol |dE| < %.0e Ha\n", ESTADO, TOL)
    println("="^96)
    @printf("%-11s | %9s %9s %7s | %8s %8s %6s | %12s\n",
            "elemento", "iter sin", "iter con", "factor", "enciende", "congela", "reini", "|dE| entre")
    @printf("%-11s | %9s %9s %7s | %8s %8s %6s | %12s\n",
            "", "DIIS", "DIIS", "", "DIIS", "DIIS", "cios", "metodos (Ha)")
    println("-"^96)
    for res in results
        i_off, i_on = res.off.r.iters, res.on.r.iters
        dE = abs(res.off.r.E_total - res.on.r.E_total)
        @printf("%-11s | %9d %9d %6.1fx | %8d %8d %6d | %12.2e\n",
                res.name, i_off, i_on, i_off / i_on,
                res.on.r.switch_iter, res.on.r.freeze_iter, res.on.r.restarts, dE)
    end
    println("-"^96)

    println("\n Energias (Ha) y contraste con el .jld2 que ya usa la tesis")
    @printf("%-11s | %18s %18s %18s | %11s %11s\n",
            "elemento", "sin DIIS", "con DIIS", "referencia", "d(sin-ref)", "d(con-ref)")
    println("-"^96)
    for res in results
        @printf("%-11s | %18.10f %18.10f %18.10f | %11.2e %11.2e\n",
                res.name, res.off.r.E_total, res.on.r.E_total, res.E_ref,
                abs(res.off.r.E_total - res.E_ref), abs(res.on.r.E_total - res.E_ref))
    end
    println("-"^96)

    println("\n Tiempo de pared (s)")
    @printf("%-11s | %12s %12s %8s\n", "elemento", "sin DIIS", "con DIIS", "factor")
    println("-"^96)
    for res in results
        @printf("%-11s | %12.1f %12.1f %7.1fx\n",
                res.name, res.off.wall, res.on.wall, res.off.wall / res.on.wall)
    end
    println("-"^96)

    bad = [r.name for r in results if abs(r.off.r.E_total - r.on.r.E_total) > 1e-8]
    if isempty(bad)
        println("\n Los dos metodos llegan a la misma energia en todos los elementos.")
    else
        println("\n AVISO: ", join(bad, ", "),
                " convergen a energias distintas. C-DIIS estaria cayendo en OTRA solucion,")
        println(" y eso es un hallazgo, no una mejora. No usar esos numeros.")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    sel = isempty(ARGS) ? ORDER : [a for a in ORDER if a in ARGS]
    isempty(sel) && error("Elementos no reconocidos: $(ARGS). Validos: $(join(ORDER, ", "))")
    results = [run_case(name; verbose = true) for name in sel]
    report(results)
    println("\n Trazas de convergencia:")
    write_traces(results; dir = @__DIR__)
end
