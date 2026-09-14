# =============================================================================
#  Codigo comun de la regeneracion de resultados de la tesis
#
#  Carga el proyecto, la configuracion (config.jl) y el manifiesto de resultados SCF
#  (benchmarks/RESULTS.toml), y reune lo que comparten las etapas: nombres de archivo,
#  espacios de trabajo, lectura y escritura de resultados_ci.toml y formato de tablas.
#  No calcula fisica por si mismo.
# =============================================================================

import Pkg
const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
Pkg.activate(REPO; io = devnull)

using AtomicSplines, JLD2, LinearAlgebra, Printf, SHA, TOML

include(joinpath(@__DIR__, "config.jl"))

const NP2           = normpath(joinpath(@__DIR__, ".."))
const MANIFIESTO    = joinpath(REPO, "benchmarks", "RESULTS.toml")
const DIR_TABLAS    = joinpath(REPO, "docs", "tablas")
const DIR_FIGURAS   = joinpath(REPO, "docs", "figures")
const RESULTADOS_CI = joinpath(NP2, "resultados_ci.toml")

# Las mismas conversiones que usan np2_ci_full.jl (AU2CM) y compute_zeta.
const HA2CM = 219474.63
const HA2EV = 27.211386245988
const ALFA  = 1.0 / 137.035999

const ELEMENTOS = ["C", "Si", "Ge", "Sn"]

# Malla de cada *_rohf.jl. `verificar_malla` la contrasta con RESULTS.toml: si un script
# cambia de malla, esto truena en vez de calcular con la geometria equivocada.
const INFO = Dict(
    "C"  => (nombre = "Carbon",    prefijo = "carbon",    etiqueta = "Carbono",
             Z = 6.0,  n_val = 2, N = 100, K = 7, gamma = 2.5),
    "Si" => (nombre = "Silicon",   prefijo = "silicon",   etiqueta = "Silicio",
             Z = 14.0, n_val = 3, N = 100, K = 7, gamma = 2.5),
    "Ge" => (nombre = "Germanium", prefijo = "germanium", etiqueta = "Germanio",
             Z = 32.0, n_val = 4, N = 300, K = 8, gamma = 3.0),
    "Sn" => (nombre = "Tin",       prefijo = "tin",       etiqueta = "Estaño",
             Z = 50.0, n_val = 5, N = 500, K = 8, gamma = 4.0),
)

orden_elemento(el) = findfirst(==(el), ELEMENTOS)

"""Valor que sigue a `nombre` en la linea de comandos, o `defecto` si no esta."""
function argumento(nombre, defecto)
    i = findfirst(==(nombre), ARGS)
    return (i === nothing || i == length(ARGS)) ? defecto : ARGS[i + 1]
end

"""Nombre del .jld2 de resultado SCF, con la nomenclatura de benchmarks/README.md."""
function archivo_scf(el::String, estado::String, alpha_d::Real = 0.0)
    suf = alpha_d > 0 ? @sprintf("_ad%.3f", alpha_d) : ""
    return "$(INFO[el].prefijo)_rohf_results_$(estado)_R$(R_MAX)$(suf).jld2"
end

ruta_np2(archivo) = joinpath(NP2, archivo)

sha256_de(ruta) = open(io -> bytes2hex(sha256(io)), ruta, "r")

corto(sha) = sha[1:8]

"""Entradas de RESULTS.toml indexadas por nombre de archivo."""
function manifiesto()
    m = TOML.parsefile(MANIFIESTO)
    return Dict(basename(r["archivo"]) => r for r in m["resultado"])
end

"""
    registrado(archivo, man) -> Bool

El archivo existe, esta en el manifiesto y su sha256 coincide. Un .jld2 sin registrar, o
cuyo sha256 cambio, no se usa: primero hay que reemitir el manifiesto y commitear.
"""
function registrado(archivo, man)
    ruta = ruta_np2(archivo)
    return isfile(ruta) && haskey(man, archivo) && sha256_de(ruta) == man[archivo]["sha256"]
end

function verificar_malla(el, archivo, man)
    haskey(man, archivo) || return nothing
    r, i = man[archivo], INFO[el]
    (r["N"] == i.N && r["K"] == i.K && r["gamma"] == i.gamma && r["R_max"] == R_MAX) ||
        error("$archivo: la malla de RESULTS.toml no coincide con INFO[\"$el\"] de comun.jl.")
    return nothing
end

"""Espacio de trabajo (geometria en cache) del elemento, con V_pol si alpha_d > 0."""
function workspace(el::String; alpha_d::Real = 0.0)
    i = INFO[el]
    return cached_init_scf_workspace(R_MAX, i.N, Val(i.K), i.Z; γ = i.gamma,
                                     alpha_d = Float64(alpha_d), r_c = R_C)
end

"""Commit corto de HEAD, con `-sucio` si hay cambios sin commitear en el codigo (.jl)."""
function commit_actual()
    h = strip(read(`git -C $REPO rev-parse --short HEAD`, String))
    cambios = strip(read(`git -C $REPO status --porcelain -- src ':(glob)benchmarks/**/*.jl'`, String))
    return isempty(cambios) ? String(h) : h * "-sucio"
end

# -----------------------------------------------------------------------------
#  resultados_ci.toml
# -----------------------------------------------------------------------------
const CABECERA_CI = """
# resultados_ci.toml - resultados del CI de valencia que alimentan las tablas de la tesis.
#
# GENERADO por benchmarks/np2_sequence/tesis/etapa2_ci.jl. No editar a mano.
#
# Cada [[ci]] es una corrida de run_full_ci (np2_ci_full.jl) sobre un .jld2 registrado en
# benchmarks/RESULTS.toml e identificado por su sha256: si el .jld2 cambia, la entrada deja
# de valer y la etapa 2 la vuelve a calcular. `m` es el numero de orbitales por cada
# l = 0..lmax. Energias en Hartree; niveles en cm^-1 desde 3P_0, en el orden
# [3P_0, 3P_1, 3P_2, 1D_2, 1S_0]. `tiempo_s` es tiempo de pared; con cache_heredada = true
# es incremental, porque ese tamano reutilizo las integrales R^k del anterior.
"""

function leer_ci(ruta = RESULTADOS_CI)
    isfile(ruta) || return Dict{String,Any}[]
    return Vector{Dict{String,Any}}(get(TOML.parsefile(ruta), "ci", Dict{String,Any}[]))
end

"""Reescribe el archivo completo, ordenado, pasando por un temporal para no truncarlo."""
function escribir_ci(entradas, ruta = RESULTADOS_CI)
    sort!(entradas, by = e -> (orden_elemento(e["elemento"]), e["estado"], e["alpha_d"], e["m"]))
    tmp = ruta * ".tmp"
    open(tmp, "w") do io
        println(io, CABECERA_CI)
        TOML.print(io, Dict("ci" => entradas); sorted = true)
    end
    mv(tmp, ruta; force = true)
    return ruta
end

"""Entrada vigente del CI para (archivo, m): la que corresponde al sha256 actual del .jld2."""
function buscar_ci(entradas, archivo, m::Int)
    ruta = ruta_np2(archivo)
    isfile(ruta) || return nothing
    sha = sha256_de(ruta)
    i = findfirst(e -> e["archivo"] == archivo && e["m"] == m && e["sha256_entrada"] == sha,
                  entradas)
    return i === nothing ? nothing : entradas[i]
end

# -----------------------------------------------------------------------------
#  Fragmentos de LaTeX
# -----------------------------------------------------------------------------
"""
    escribir_tabla(nombre, cuerpo, procedencia; dir)

Escribe un fragmento `tabular` para `\\input{tablas/<nombre>}`. La leyenda y la etiqueta
se quedan en el capitulo: aqui solo van numeros y su procedencia. La emision es
reproducible, sin fecha ni commit: si los numeros no cambian, el archivo tampoco.
"""
function escribir_tabla(nombre, cuerpo, procedencia; dir = DIR_TABLAS)
    mkpath(dir)
    ruta = joinpath(dir, nombre)
    open(ruta, "w") do io
        println(io, "% GENERADO por benchmarks/np2_sequence/tesis/etapa3_tablas.jl; no editar a mano.")
        println(io, "% Procedencia:")
        foreach(l -> println(io, "%   ", l), procedencia)
        print(io, cuerpo)
    end
    println("  tabla -> ", startswith(ruta, REPO) ? relpath(ruta, REPO) : ruta)
    return ruta
end

# -----------------------------------------------------------------------------
#  Seleccion de alpha_d y curvas de convergencia (las usan las etapas 3 y 4)
# -----------------------------------------------------------------------------
error_3P(niv, el) = (niv[2:3] .- NIST_NIVELES[el][2:3]) ./ NIST_NIVELES[el][2:3]
costo_alpha(niv, el) = CRITERIO_ALPHA == :solo_3P2 ? abs(error_3P(niv, el)[2]) :
                                                     sqrt(sum(abs2, error_3P(niv, el)) / 2)

"""
    alpha_elegido(el, ci) -> (alpha_d, entrada) o nothing

alpha_d del barrido con menor costo segun CRITERIO_ALPHA. `ci(el, estado, alpha_d)` devuelve la
entrada de resultados_ci.toml que corresponde, o nothing.
"""
function alpha_elegido(el, ci)
    cands = [(a, ci(el, ESTADO, a)) for a in get(BARRIDO_ALPHA, el, Float64[])]
    filter!(c -> c[2] !== nothing, cands)
    isempty(cands) && return nothing
    return argmin(c -> costo_alpha(c[2]["niveles_cm"], el), cands)
end

"""Entradas vigentes de la curva de convergencia (ESTADO, alpha_d = 0) del elemento, por m."""
function curva_convergencia(ci, man, el)
    archivo = archivo_scf(el, ESTADO, 0.0)
    registrado(archivo, man) || return Dict{String,Any}[]
    sha = sha256_de(ruta_np2(archivo))
    return sort([e for e in ci if e["archivo"] == archivo && e["sha256_entrada"] == sha]; by = e -> e["m"])
end
