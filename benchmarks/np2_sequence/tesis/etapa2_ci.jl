# =============================================================================
#  Etapa 2: CI de valencia (la parte larga)
#
#  Corre run_full_ci (np2_ci_full.jl) sobre los .jld2 registrados y guarda cada resultado
#  en benchmarks/np2_sequence/resultados_ci.toml en cuanto termina. Si la corrida se
#  interrumpe no se pierde lo ya calculado, y al relanzarla se salta lo que ya esta.
#
#  Grupos:
#    convergencia  ESTADO, alpha_d = 0, m en TAMANOS_CONVERGENCIA. El ultimo tamano es el
#                  de produccion: el HF+CI de las tablas.
#    sensibilidad  ESTADO_SENSIBILIDAD, alpha_d = 0, m = M_PRODUCCION.
#    barrido       ESTADO, cada alpha_d de BARRIDO_ALPHA, m = M_PRODUCCION.
#
#  Uso, desde la raiz del repo (del orden de 2-3 h en total):
#    julia --project=. benchmarks/np2_sequence/tesis/etapa2_ci.jl --lista
#    julia --project=. benchmarks/np2_sequence/tesis/etapa2_ci.jl 2>&1 | tee etapa2.log
#    ... --solo Ge,Sn --grupo barrido           filtros combinables
#    ... --prueba --salida /ruta/prueba.toml    espacio minimo, para ver que corre
# =============================================================================

include(joinpath(@__DIR__, "comun.jl"))
include(joinpath(NP2, "np2_ci_full.jl"))

# Tiempo de pared a m = 20, lmax = 3 (HALLAZGOS-2026-09-06 §4); solo para --lista.
const SEGUNDOS_M20 = Dict("C" => 92.0, "Si" => 83.0, "Ge" => 390.0, "Sn" => 609.0)

const PRUEBA  = "--prueba" in ARGS
const TAMANOS = PRUEBA ? [2, 3] : TAMANOS_CONVERGENCIA
const M_PROD  = last(TAMANOS)
const SALIDA  = argumento("--salida", PRUEBA ? joinpath(tempdir(), "resultados_ci_prueba.toml") :
                                               RESULTADOS_CI)
const SOLO    = split(argumento("--solo", join(ELEMENTOS, ",")), ",")
const GRUPOS  = split(argumento("--grupo", "convergencia,sensibilidad,barrido"), ",")

"""Casos a correr: (grupo, elemento, estado, alpha_d, tamanos)."""
function casos()
    cs = Tuple{String,String,String,Float64,Vector{Int}}[]
    for el in ELEMENTOS
        el in SOLO || continue
        "convergencia" in GRUPOS && push!(cs, ("convergencia", el, ESTADO, 0.0, TAMANOS))
        ("sensibilidad" in GRUPOS && ESTADO_SENSIBILIDAD != ESTADO) &&
            push!(cs, ("sensibilidad", el, ESTADO_SENSIBILIDAD, 0.0, [M_PROD]))
        if "barrido" in GRUPOS
            for a in get(BARRIDO_ALPHA, el, Float64[])
                push!(cs, ("barrido", el, ESTADO, a, [M_PROD]))
            end
        end
    end
    return cs
end

function entrada(el, estado, archivo, sha, alpha_d, m, r, t, verificado, heredada, commit)
    return Dict{String,Any}(
        "elemento" => el, "estado" => estado, "archivo" => archivo, "sha256_entrada" => sha,
        "alpha_d" => alpha_d, "r_c" => R_C, "m" => m, "lmax" => LMAX_CI,
        "n_orb" => r.n_orb, "n_csf" => collect(r.n_csf), "n_rk" => r.n_rk,
        "zeta_Ha" => r.zeta, "F2_Ha" => r.F2, "E_ref_Ha" => r.E_3P - r.E_corr,
        "E_3P_Ha" => r.E_3P, "E_1D_Ha" => r.E_1D, "E_1S_Ha" => r.E_1S,
        "E_corr_Ha" => r.E_corr, "pesos_np2" => collect(r.C), "raices" => collect(r.roots),
        "g_eff_3P2" => r.g_eff, "niveles_cm" => r.levels,
        # autoverificacion = true: T1-T4 corrieron sobre este tamano. false: compartio workspace
        # con un tamano menor del mismo archivo que ya las paso en esta corrida.
        "autoverificacion" => verificado,
        # cache_heredada = true: tiempo_s es incremental, porque reutilizo las R^k del tamano
        # anterior. false: primer tamano, o corrida reanudada, y el tiempo incluye todas.
        "cache_heredada" => heredada,
        "tiempo_s" => round(t; digits = 1), "commit" => commit)
end

function main_etapa2()
    man = manifiesto()
    entradas = leer_ci(SALIDA)
    commit = commit_actual()
    (endswith(commit, "-sucio") && !PRUEBA) &&
        @warn "Hay cambios de codigo sin commitear: las entradas quedaran marcadas $commit."
    cs = casos()

    println("\nCasos$(PRUEBA ? " (PRUEBA)" : "") -> $(SALIDA)")
    pendiente = 0.0
    for (grupo, el, estado, a, ms) in cs
        archivo = archivo_scf(el, estado, a)
        faltan = [m for m in ms if buscar_ci(entradas, archivo, m) === nothing]
        ok = registrado(archivo, man)
        seg = isempty(faltan) ? 0.0 : SEGUNDOS_M20[el] * (length(faltan) > 1 ? 1.15 : 1.0)
        ok && (pendiente += seg)
        estado_txt = !ok ? "SIN REGISTRAR (etapa 1 y reemitir el manifiesto)" :
                     isempty(faltan) ? "listo" : @sprintf("~%.0f min", seg / 60)
        @printf("  %-13s %-46s m = %-15s %s\n", grupo, archivo, join(faltan, ","), estado_txt)
    end
    PRUEBA || @printf("Pendiente: del orden de %.0f min.\n", pendiente / 60)
    "--lista" in ARGS && return

    for (grupo, el, estado, a, ms) in cs
        archivo = archivo_scf(el, estado, a)
        if !registrado(archivo, man)
            @warn "Se omite $archivo: no esta registrado en benchmarks/RESULTS.toml."
            continue
        end
        verificar_malla(el, archivo, man)
        # run_full_ci toma alpha_d del .jld2 con get(data, "alpha_d", 0.0): un archivo con V_pol
        # que no guardara la clave se calcularia en silencio sobre la geometria sin V_pol.
        a_jld2 = jldopen(f -> haskey(f, "alpha_d") ? f["alpha_d"] : 0.0, ruta_np2(archivo), "r")
        a_jld2 == man[archivo]["alpha_d"] ||
            error("$archivo: alpha_d = $a_jld2 en el .jld2 y $(man[archivo]["alpha_d"]) en RESULTS.toml.")
        sha = sha256_de(ruta_np2(archivo))
        ws = nothing
        verificado = false
        for m in ms
            buscar_ci(entradas, archivo, m) === nothing || continue
            @printf("\n##### %s | %s | m = %d, lmax = %d #####\n", grupo, archivo, m, LMAX_CI)
            heredada = ws !== nothing
            t0 = time()
            r = run_full_ci(INFO[el].nombre, archivo; n_per_l = espacio_activo(m),
                            selftest = !verificado, verbose = false, ws = ws)
            t = time() - t0
            filter!(e -> !(e["archivo"] == archivo && e["m"] == m), entradas)
            push!(entradas, entrada(el, estado, archivo, sha, a, m, r, t, !verificado, heredada, commit))
            escribir_ci(entradas, SALIDA)
            verificado = true
            ws = r.ws        # cache de R^k compartida con el siguiente tamano del mismo archivo
            @printf("  -> E_corr = %.6e Ha | 3P_2 = %.2f cm^-1 | %.1f s | guardado\n",
                    r.E_corr, r.levels[3], t)
        end
    end
    println("\nEtapa 2 terminada. Siguen etapa3_tablas.jl y etapa4_figuras.jl.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_etapa2()
end
