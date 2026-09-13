#!/usr/bin/env julia
#
#     julia benchmarks/verificar_resultados.jl            # verifica
#     julia --project=. benchmarks/verificar_resultados.jl --emitir > benchmarks/RESULTS.toml
#
# Verifica que los .jld2 de resultado SCF de benchmarks/ sigan siendo los que
# produjeron las tablas del capitulo de resultados, comparando el sha256 de cada
# archivo contra benchmarks/RESULTS.toml.
#
# Por que existe: el 2026-09-06 se descubrio que el signo de coeff_k2 del termino
# 3P se habia invertido en carbon_rohf.jl y silicon_rohf.jl *despues* de generar
# las tablas del capitulo 6. Los .jld2 no estaban trackeados, asi que no hubo
# forma de fechar la regresion con git y hubo que reconstruirla comparando a mano
# contra Froese Fischer. Con el manifiesto, un cambio en cualquier resultado
# aparece como un sha256 distinto y el diff del .toml dice exactamente que
# energia se movio y en cuanto.
#
# El modo de verificacion NO carga el proyecto: solo usa SHA y TOML, que son de
# la biblioteca estandar, para que pueda correrse en cualquier parte (un hook de
# git, CI) sin instanciar el entorno. El modo --emitir si lo carga, porque lee
# E_total de los .jld2 con JLD2.
#
# Los metadatos (parametros fisicos, script generador, notas) viven en el
# CATALOGO de abajo, no en el .toml: asi regenerar el manifiesto no borra lo que
# se sabe de cada archivo. Al registrar un resultado nuevo, agregalo aqui.

using SHA, TOML, Printf

const BENCH      = @__DIR__
const RAIZ       = normpath(joinpath(BENCH, ".."))
const MANIFIESTO = joinpath(BENCH, "RESULTS.toml")

# (Z, R_max, N, K, gamma, alpha_d, r_c) es la clave de la geometria; `estado`
# distingue el promedio de configuracion del termino 3P especifico.
const CATALOGO = [
    (archivo = "np2_sequence/carbon_rohf_results_av_R30.0.jld2",
     elemento = "C", Z = 6.0, estado = "av", R_max = 30.0, N = 100, K = 7,
     gamma = 2.5, alpha_d = 0.0, r_c = 1.0,
     script = "np2_sequence/carbon_rohf.jl",
     nota = "Promedio de configuracion. coeff_k2 = +2/25, correcto desde siempre."),

    (archivo = "np2_sequence/carbon_rohf_results_3P_R30.0.jld2",
     elemento = "C", Z = 6.0, estado = "3P", R_max = 30.0, N = 100, K = 7,
     gamma = 2.5, alpha_d = 0.0, r_c = 1.0,
     script = "np2_sequence/carbon_rohf.jl",
     nota = "Limite HF: coincide con Froese Fischer y con resultados.tex:25. " *
            "Reemplazo al archivo contaminado por el signo de f2 (ver HALLAZGOS-2026-09-06). " *
            "PENDIENTE: se escribio a mano, le faltan las claves active_s y active_p; " *
            "regenerar con el script actual (segundos) para completar el esquema."),

    (archivo = "np2_sequence/silicon_rohf_results_av_R30.0.jld2",
     elemento = "Si", Z = 14.0, estado = "av", R_max = 30.0, N = 100, K = 7,
     gamma = 2.5, alpha_d = 0.0, r_c = 1.0,
     script = "np2_sequence/silicon_rohf.jl",
     nota = "Promedio de configuracion. coeff_k2 = +2/25, correcto desde siempre."),

    (archivo = "np2_sequence/silicon_rohf_results_3P_R30.0.jld2",
     elemento = "Si", Z = 14.0, estado = "3P", R_max = 30.0, N = 100, K = 7,
     gamma = 2.5, alpha_d = 0.0, r_c = 1.0,
     script = "np2_sequence/silicon_rohf.jl",
     nota = "Limite HF: coincide con Froese Fischer y con resultados.tex:30. " *
            "Reemplazo al archivo contaminado por el signo de f2 (ver HALLAZGOS-2026-09-06). " *
            "PENDIENTE: le faltan las claves active_s y active_p, igual que en carbono."),

    (archivo = "np2_sequence/germanium_rohf_results_av_R30.0.jld2",
     elemento = "Ge", Z = 32.0, estado = "av", R_max = 30.0, N = 300, K = 8,
     gamma = 3.0, alpha_d = 0.0, r_c = 1.0,
     script = "np2_sequence/germanium_rohf.jl",
     nota = "Promedio de configuracion."),

    (archivo = "np2_sequence/germanium_rohf_results_3P_R30.0.jld2",
     elemento = "Ge", Z = 32.0, estado = "3P", R_max = 30.0, N = 300, K = 8,
     gamma = 3.0, alpha_d = 0.0, r_c = 1.0,
     script = "np2_sequence/germanium_rohf.jl",
     nota = "Nunca afectado por el signo de f2: germanio siempre uso +5/25."),

    (archivo = "np2_sequence/germanium_rohf_results_av_R30.0_ad0.500.jld2",
     elemento = "Ge", Z = 32.0, estado = "av", R_max = 30.0, N = 300, K = 8,
     gamma = 3.0, alpha_d = 0.5, r_c = 1.0,
     script = "np2_sequence/germanium_rohf_vpol.jl",
     nota = "Con potencial de polarizacion del core. Es la columna CI+V_pol del " *
            "Cuadro tab:niveles_energia_pesados (resultados.tex:248). germanium_rohf.jl con " *
            "alpha_d = 0.5 y C-DIIS lo reproduce a 1.4e-8 Ha."),

    (archivo = "np2_sequence/tin_rohf_results_3P_R30.0.jld2",
     elemento = "Sn", Z = 50.0, estado = "3P", R_max = 30.0, N = 500, K = 8,
     gamma = 4.0, alpha_d = 0.0, r_c = 1.0,
     script = "np2_sequence/tin_rohf.jl",
     nota = "Limite HF del estanio. Nunca afectado por el signo de f2."),

    (archivo = "np2_sequence/tin_rohf_results_av_R30.0.jld2",
     elemento = "Sn", Z = 50.0, estado = "av", R_max = 30.0, N = 500, K = 8,
     gamma = 4.0, alpha_d = 0.0, r_c = 1.0,
     script = "np2_sequence/tin_rohf.jl",
     nota = "Promedio de configuracion del estanio, el que faltaba en toda la secuencia. " *
            "E(3P) - E(av) = -0.017531 Ha contra la prediccion -0.12*F2 = -0.017668, " *
            "un 0.8% de diferencia por relajacion orbital, igual que en germanio. " *
            "Es la columna vacia del estanio en tab:niveles_energia_pesados."),

    (archivo = "np2_sequence/tin_rohf_results_3P_R30.0_ad1.000.jld2",
     elemento = "Sn", Z = 50.0, estado = "3P", R_max = 30.0, N = 500, K = 8,
     gamma = 4.0, alpha_d = 1.0, r_c = 1.0,
     script = "np2_sequence/tin_rohf_vpol.jl",
     nota = "Generado por tin_rohf_vpol.jl (rama 3P con +5/25, correcta). tin_rohf.jl con alpha_d = 1 " *
            "y C-DIIS lo reproduce a 1.4e-8 Ha. Es el punto alpha_d = 1 del barrido 3P de " *
            "tesis/config.jl y la entrada de tin_fine_structure.jl:252."),

    (archivo = "np2_sequence/tin_rohf_results_av_R30.0_ad1.000.jld2",
     elemento = "Sn", Z = 50.0, estado = "av", R_max = 30.0, N = 500, K = 8,
     gamma = 4.0, alpha_d = 1.0, r_c = 1.0,
     script = "np2_sequence/tin_rohf_vpol.jl",
     nota = "Barrido de alpha_d del estanio. Este valor quedo fuera del Cuadro " *
            "tab:estanio_ci, que arranca en alpha_d = 2."),

    (archivo = "np2_sequence/tin_rohf_results_av_R30.0_ad2.000.jld2",
     elemento = "Sn", Z = 50.0, estado = "av", R_max = 30.0, N = 500, K = 8,
     gamma = 4.0, alpha_d = 2.0, r_c = 1.0,
     script = "np2_sequence/tin_rohf_vpol.jl",
     nota = "Renglon alpha_d = 2.0 del Cuadro tab:estanio_ci (resultados.tex:281)."),

    (archivo = "np2_sequence/tin_rohf_results_av_R30.0_ad3.000.jld2",
     elemento = "Sn", Z = 50.0, estado = "av", R_max = 30.0, N = 500, K = 8,
     gamma = 4.0, alpha_d = 3.0, r_c = 1.0,
     script = "np2_sequence/tin_rohf_vpol.jl",
     nota = "Renglon alpha_d = 3.0 del Cuadro tab:estanio_ci (resultados.tex:281)."),

    (archivo = "np2_sequence/tin_rohf_results_av_R30.0_ad4.000.jld2",
     elemento = "Sn", Z = 50.0, estado = "av", R_max = 30.0, N = 500, K = 8,
     gamma = 4.0, alpha_d = 4.0, r_c = 1.0,
     script = "np2_sequence/tin_rohf_vpol.jl",
     nota = "Renglon alpha_d = 4.0 del Cuadro tab:estanio_ci (resultados.tex:281)."),

    (archivo = "np2_sequence/tin_rohf_results_av_R30.0_ad5.000.jld2",
     elemento = "Sn", Z = 50.0, estado = "av", R_max = 30.0, N = 500, K = 8,
     gamma = 4.0, alpha_d = 5.0, r_c = 1.0,
     script = "np2_sequence/tin_rohf_vpol.jl",
     nota = "Renglon alpha_d = 5.0 del Cuadro tab:estanio_ci (resultados.tex:281)."),

    (archivo = "np2_sequence/tin_rohf_results_av_R30.0_ad6.000.jld2",
     elemento = "Sn", Z = 50.0, estado = "av", R_max = 30.0, N = 500, K = 8,
     gamma = 4.0, alpha_d = 6.0, r_c = 1.0,
     script = "np2_sequence/tin_rohf_vpol.jl",
     nota = "Renglon alpha_d = 6.0 del Cuadro tab:estanio_ci (resultados.tex:281)."),

    # Planeados por benchmarks/np2_sequence/tesis/config.jl. Hasta que la etapa 1 los genere,
    # --emitir los omite con un aviso.
    (archivo = "np2_sequence/germanium_rohf_results_3P_R30.0_ad0.250.jld2",
     elemento = "Ge", Z = 32.0, estado = "3P", R_max = 30.0, N = 300, K = 8,
     gamma = 3.0, alpha_d = 0.25, r_c = 1.0,
     script = "np2_sequence/germanium_rohf.jl",
     nota = "Barrido de V_pol sobre orbitales 3P de tesis/config.jl. germanium_rohf.jl con alpha_d y C-DIIS."),

    (archivo = "np2_sequence/germanium_rohf_results_3P_R30.0_ad0.500.jld2",
     elemento = "Ge", Z = 32.0, estado = "3P", R_max = 30.0, N = 300, K = 8,
     gamma = 3.0, alpha_d = 0.5, r_c = 1.0,
     script = "np2_sequence/germanium_rohf.jl",
     nota = "Barrido de V_pol sobre orbitales 3P de tesis/config.jl. germanium_rohf.jl con alpha_d y C-DIIS."),

    (archivo = "np2_sequence/germanium_rohf_results_3P_R30.0_ad0.750.jld2",
     elemento = "Ge", Z = 32.0, estado = "3P", R_max = 30.0, N = 300, K = 8,
     gamma = 3.0, alpha_d = 0.75, r_c = 1.0,
     script = "np2_sequence/germanium_rohf.jl",
     nota = "Barrido de V_pol sobre orbitales 3P de tesis/config.jl. germanium_rohf.jl con alpha_d y C-DIIS."),

    (archivo = "np2_sequence/germanium_rohf_results_3P_R30.0_ad1.000.jld2",
     elemento = "Ge", Z = 32.0, estado = "3P", R_max = 30.0, N = 300, K = 8,
     gamma = 3.0, alpha_d = 1.0, r_c = 1.0,
     script = "np2_sequence/germanium_rohf.jl",
     nota = "Barrido de V_pol sobre orbitales 3P de tesis/config.jl. germanium_rohf.jl con alpha_d y C-DIIS."),

    (archivo = "np2_sequence/tin_rohf_results_3P_R30.0_ad2.000.jld2",
     elemento = "Sn", Z = 50.0, estado = "3P", R_max = 30.0, N = 500, K = 8,
     gamma = 4.0, alpha_d = 2.0, r_c = 1.0,
     script = "np2_sequence/tin_rohf.jl",
     nota = "Barrido de V_pol sobre orbitales 3P de tesis/config.jl. tin_rohf.jl con alpha_d y C-DIIS."),

    (archivo = "np2_sequence/tin_rohf_results_3P_R30.0_ad3.000.jld2",
     elemento = "Sn", Z = 50.0, estado = "3P", R_max = 30.0, N = 500, K = 8,
     gamma = 4.0, alpha_d = 3.0, r_c = 1.0,
     script = "np2_sequence/tin_rohf.jl",
     nota = "Barrido de V_pol sobre orbitales 3P de tesis/config.jl. tin_rohf.jl con alpha_d y C-DIIS."),

    (archivo = "np2_sequence/tin_rohf_results_3P_R30.0_ad4.000.jld2",
     elemento = "Sn", Z = 50.0, estado = "3P", R_max = 30.0, N = 500, K = 8,
     gamma = 4.0, alpha_d = 4.0, r_c = 1.0,
     script = "np2_sequence/tin_rohf.jl",
     nota = "Barrido de V_pol sobre orbitales 3P de tesis/config.jl. tin_rohf.jl con alpha_d y C-DIIS."),

    (archivo = "np2_sequence/tin_rohf_results_3P_R30.0_ad5.000.jld2",
     elemento = "Sn", Z = 50.0, estado = "3P", R_max = 30.0, N = 500, K = 8,
     gamma = 4.0, alpha_d = 5.0, r_c = 1.0,
     script = "np2_sequence/tin_rohf.jl",
     nota = "Barrido de V_pol sobre orbitales 3P de tesis/config.jl. tin_rohf.jl con alpha_d y C-DIIS."),

    (archivo = "np2_sequence/tin_rohf_results_3P_R30.0_ad6.000.jld2",
     elemento = "Sn", Z = 50.0, estado = "3P", R_max = 30.0, N = 500, K = 8,
     gamma = 4.0, alpha_d = 6.0, r_c = 1.0,
     script = "np2_sequence/tin_rohf.jl",
     nota = "Barrido de V_pol sobre orbitales 3P de tesis/config.jl. tin_rohf.jl con alpha_d y C-DIIS."),
]

sha_de(ruta) = open(ruta, "r") do io; bytes2hex(sha256(io)); end

function git(cmd::Cmd)
    try
        return strip(read(setenv(cmd, dir = RAIZ), String))
    catch
        return ""
    end
end

# --------------------------------------------------------------------------
# Verificacion
# --------------------------------------------------------------------------
function verificar()
    isfile(MANIFIESTO) || (println("No existe $MANIFIESTO."); return 1)
    manif = TOML.parsefile(MANIFIESTO)
    registros = get(manif, "resultado", Dict{String,Any}[])
    fallas = String[]
    registrados = Set{String}()

    for r in registros
        archivo = r["archivo"]
        push!(registrados, archivo)
        ruta = joinpath(BENCH, archivo)
        if !isfile(ruta)
            push!(fallas, "FALTA      $archivo")
            continue
        end
        h = sha_de(ruta)
        if h == r["sha256"]
            @printf("ok         %-52s E = %s\n", archivo, string(r["E_total"]))
        else
            push!(fallas, "CAMBIO     $archivo\n" *
                          "           registrado: $(r["sha256"])\n" *
                          "           en disco:   $h\n" *
                          "           E_total registrada: $(r["E_total"]) ($(r["elemento"]) $(r["estado"]))")
        end
    end

    # Archivos de resultado que nadie registro.
    for (raiz, _, nombres) in walkdir(BENCH), n in nombres
        occursin("rohf_results", n) && endswith(n, ".jld2") || continue
        rel = relpath(joinpath(raiz, n), BENCH)
        rel in registrados || push!(fallas, "SIN REGISTRAR  $rel")
    end

    println()
    if isempty(fallas)
        println("Los $(length(registros)) resultados coinciden con el manifiesto.")
        return 0
    end
    println("$(length(fallas)) problema(s):\n")
    for f in fallas; println(f); end
    println("""

    Un CAMBIO significa que el .jld2 ya no es el que produjo las tablas. Si fue a
    proposito (regeneraste el SCF), reemitir el manifiesto y commitear el diff:

        julia --project=. benchmarks/verificar_resultados.jl --emitir > benchmarks/RESULTS.toml

    Si no fue a proposito, `git checkout` sobre el archivo lo devuelve.
    """)
    return 1
end

# --------------------------------------------------------------------------
# Emision del manifiesto
# --------------------------------------------------------------------------
function emitir()
    println("""
    # RESULTS.toml - manifiesto de los resultados SCF que alimentan la tesis.
    #
    # GENERADO por benchmarks/verificar_resultados.jl --emitir. No editar a mano:
    # los metadatos viven en el CATALOGO de ese script.
    #
    # Cada entrada fija el sha256 del .jld2, su energia total y los parametros
    # fisicos con los que se genero. Un `git diff` sobre este archivo dice que
    # resultado se movio y en cuanto, aunque el binario no diga nada legible.
    #
    # La emision es reproducible: no lleva fecha ni sello de commit, asi que
    #
    #     julia --project=. benchmarks/verificar_resultados.jl --emitir | diff - benchmarks/RESULTS.toml
    #
    # sale vacio mientras nada haya cambiado de verdad. Cuando se reemitio y desde
    # que commit lo dice el historial de git de este mismo archivo.
    """)

    for e in CATALOGO
        ruta = joinpath(BENCH, e.archivo)
        if !isfile(ruta)
            println(stderr, "aviso: no existe $(e.archivo), se omite")
            continue
        end
        E = jldopen(ruta, "r") do f; f["E_total"]; end
        nspl = jldopen(ruta, "r") do f; haskey(f, "num_splines") ? f["num_splines"] : -1; end
        commit_codigo = git(`git log -1 --format=%h -- $(joinpath("benchmarks", e.script))`)
        sucio = !isempty(git(`git diff --name-only -- $(joinpath("benchmarks", e.script))`))

        println("[[resultado]]")
        println("archivo       = ", repr(e.archivo))
        println("elemento      = ", repr(e.elemento))
        println("estado        = ", repr(e.estado))
        println("sha256        = ", repr(sha_de(ruta)))
        println("bytes         = ", filesize(ruta))
        println("E_total       = ", E)
        println("num_splines   = ", nspl)
        println("Z             = ", e.Z)
        println("R_max         = ", e.R_max)
        println("N             = ", e.N)
        println("K             = ", e.K)
        println("gamma         = ", e.gamma)
        println("alpha_d       = ", e.alpha_d)
        println("r_c           = ", e.r_c)
        println("script        = ", repr(e.script))
        println("commit_codigo = ", repr(commit_codigo))
        sucio && println("codigo_sucio  = true  # el script tiene cambios sin commitear")
        println("nota          = ", repr(e.nota))
        println()
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    if "--emitir" in ARGS
        import Pkg
        Pkg.activate(RAIZ; io = devnull)
        @eval using JLD2
        emitir()
    else
        exit(verificar())
    end
end
