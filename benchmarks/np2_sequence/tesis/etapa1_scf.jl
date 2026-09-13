# =============================================================================
#  Etapa 1: que existan y esten registrados los .jld2 que pide la configuracion
#
#  Lista los resultados SCF que necesita config.jl y el estado de cada uno:
#    registrado     existe y su sha256 coincide con benchmarks/RESULTS.toml
#    sin_registrar  existe pero no esta en el manifiesto, o cambio: hay que reemitirlo
#    falta          hay que correr el SCF
#
#  Uso, desde la raiz del repo:
#    julia --project=. benchmarks/np2_sequence/tesis/etapa1_scf.jl                 # solo lista
#    julia --project=. benchmarks/np2_sequence/tesis/etapa1_scf.jl --prueba-vpol   # ~6 min
#    julia --project=. benchmarks/np2_sequence/tesis/etapa1_scf.jl --correr 2>&1 | tee etapa1.log
#    ... --correr --sin-diis     todo sin C-DIIS, con level shift permanente (~2.5 h)
#
#  --prueba-vpol reproduce, sin guardar, dos archivos registrados que generaron los
#  *_rohf_vpol.jl legados: germanium_rohf_results_av_R30.0_ad0.500.jld2 con germanium_rohf.jl y
#  tin_rohf_results_3P_R30.0_ad1.000.jld2 con tin_rohf.jl, los dos con alpha_d y C-DIIS. Es la
#  prueba de que la plomeria de V_pol de los scripts base equivale a la de los legados.
#
#  Con --correr, un punto que no converja con C-DIIS en 300 iteraciones se repite sin C-DIIS.
#
#  Los SCF nuevos usan germanium_rohf.jl y tin_rohf.jl con alpha_d, NO los *_rohf_vpol.jl.
# =============================================================================

include(joinpath(@__DIR__, "comun.jl"))
for s in ("carbon_rohf.jl", "silicon_rohf.jl", "germanium_rohf.jl", "tin_rohf.jl")
    include(joinpath(NP2, s))
end

# Minutos de pared por SCF con C-DIIS (HALLAZGOS-2026-09-06 §6); sin C-DIIS, de 2 a 5 veces mas.
const MINUTOS_SCF = Dict("C" => 0.1, "Si" => 0.1, "Ge" => 2.5, "Sn" => 6.0)

function requeridos()
    req = Tuple{String,String,Float64}[]
    for el in ELEMENTOS, estado in unique([ESTADO, ESTADO_SENSIBILIDAD])
        push!(req, (el, estado, 0.0))
    end
    for el in ELEMENTOS, a in get(BARRIDO_ALPHA, el, Float64[])
        push!(req, (el, ESTADO, a))
    end
    return req
end

function situacion(archivo, man)
    isfile(ruta_np2(archivo)) || return "falta"
    return registrado(archivo, man) ? "registrado" : "sin_registrar"
end

function correr_scf(el, estado, alpha_d; use_diis::Bool, max_iter::Int = 1000)
    kw = (estado = estado, use_diis = use_diis, save = true, verbose = true, max_iter = max_iter)
    if el in ("C", "Si")
        alpha_d == 0 || error("$(INFO[el].prefijo)_rohf.jl no acepta alpha_d.")
        return el == "C" ? solve_carbon_rohf(R_MAX; kw...) : solve_silicon_rohf(R_MAX; kw...)
    elseif el == "Ge"
        return solve_germanium_rohf(R_MAX; kw..., alpha_d = alpha_d, r_c = R_C)
    else
        return solve_tin_rohf(R_MAX; kw..., alpha_d = alpha_d, r_c = R_C)
    end
end

function main_etapa1()
    man = manifiesto()
    catalogo = read(joinpath(REPO, "benchmarks", "verificar_resultados.jl"), String)

    println("\nConfiguracion: ESTADO = $ESTADO, sensibilidad = $ESTADO_SENSIBILIDAD, r_c = $R_C")
    faltan = Tuple{String,String,Float64}[]
    for (el, estado, a) in requeridos()
        archivo = archivo_scf(el, estado, a)
        st = situacion(archivo, man)
        st == "falta" && push!(faltan, (el, estado, a))
        nota = occursin(archivo, catalogo) ? "" : "   <- agregarlo al CATALOGO de verificar_resultados.jl"
        @printf("  %-13s %s%s\n", st, archivo, nota)
    end
    @printf("\n%d SCF por correr, del orden de %.0f min con C-DIIS.\n", length(faltan),
            sum((MINUTOS_SCF[el] for (el, _, _) in faltan); init = 0.0))

    if "--prueba-vpol" in ARGS
        # Cada script base con alpha_d contra un archivo registrado que genero su legado
        # *_rohf_vpol.jl. Medido el 2026-09-13: 1.35e-8 Ha (Ge, 72 iter.) y 1.36e-8 Ha (Sn, 41).
        for archivo in ("germanium_rohf_results_av_R30.0_ad0.500.jld2",
                        "tin_rohf_results_3P_R30.0_ad1.000.jld2")
            reg = man[archivo]
            el = reg["elemento"]
            kw = (estado = reg["estado"], use_diis = true, save = false, verbose = false,
                  alpha_d = reg["alpha_d"], r_c = reg["r_c"])
            t0 = time()
            r = el == "Ge" ? solve_germanium_rohf(R_MAX; kw...) : solve_tin_rohf(R_MAX; kw...)
            d = r.E_total - reg["E_total"]
            @printf("\nPrueba V_pol, %s: E = %.10f | registrada = %.10f | d = %.2e Ha | %d iter | %.1f min\n",
                    archivo, r.E_total, reg["E_total"], d, r.iters, (time() - t0) / 60)
            (r.converged && abs(d) < 1e-7) ||
                error("$(INFO[el].prefijo)_rohf.jl con alpha_d NO reproduce $archivo.")
        end
        println("OK: con alpha_d, los scripts base reproducen lo que generaron los legados.")
    end

    if "--correr" in ARGS && !isempty(faltan)
        con_diis = !("--sin-diis" in ARGS)
        for (i, (el, estado, a)) in enumerate(faltan)
            archivo = archivo_scf(el, estado, a)
            @printf("\n##### [%d/%d] %s #####\n", i, length(faltan), archivo)
            t0 = time()
            # Con C-DIIS estos SCF convergen en 40-90 iteraciones. El tope de 300 acota lo que se
            # pierde si alguno no converge, y ese punto se repite sin C-DIIS con el tope de siempre.
            modo = con_diis ? "C-DIIS" : "sin C-DIIS"
            r = correr_scf(el, estado, a; use_diis = con_diis, max_iter = con_diis ? 300 : 1000)
            if con_diis && !r.converged
                @warn "$archivo: no convergio con C-DIIS en 300 iteraciones; se repite sin C-DIIS."
                modo = "sin C-DIIS, tras fallar con C-DIIS"
                r = correr_scf(el, estado, a; use_diis = false, max_iter = 1000)
            end
            (r.converged && isfile(ruta_np2(archivo))) ||
                error("$archivo: el SCF no convergio o no escribio el archivo.")
            @printf("  -> E = %.10f Ha | %d iteraciones | %s | %.1f min\n",
                    r.E_total, r.iters, modo, (time() - t0) / 60)
        end
        println("""

        Listo. Para registrar los archivos nuevos, desde la raiz del repo:
          1. Si alguno no esta en el CATALOGO de benchmarks/verificar_resultados.jl, agregarlo.
          2. julia --project=. benchmarks/verificar_resultados.jl --emitir > benchmarks/RESULTS.toml
          3. julia benchmarks/verificar_resultados.jl
          4. Commitear los .jld2 nuevos y el manifiesto.
        La etapa 2 no usa archivos sin registrar.""")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_etapa1()
end
