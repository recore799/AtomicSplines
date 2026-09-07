using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using AtomicSplines
using JLD2
using Printf

include("np2_toy_ci.jl")

# -----------------------------------------------------------------------------
# CONFIGURACION -- rellenar a mano y a conciencia.
#
# La version anterior de este script llamaba run_toy_ci(Float64(Z), alpha_d, r_c),
# una firma que ya no existe: fallaba con MethodError. Es decir, las tablas del
# capitulo 6 no se podian regenerar. Aqui las rutas son explicitas para que la
# procedencia de cada numero quede escrita en el script y no en la memoria de nadie.
#
#   file_hf   : .jld2 del SCF sin V_pol
#   file_vpol : .jld2 del SCF con V_pol (nothing => la columna repite la de HF+CI)
#
# zeta ya NO se pasa como literal: run_toy_ci lo calcula de R_grid/V_eff/P_np del
# propio archivo, que es lo unico que hace que un barrido de alpha_d tenga sentido.
# -----------------------------------------------------------------------------
const CASES = [
    (element = "Carbon",    label = "Carbono (C I)",  Z =  6,
     file_hf = "carbon_rohf_results_3P_R30.0.jld2",     file_vpol = nothing),
    (element = "Silicon",   label = "Silicio (Si I)", Z = 14,
     file_hf = "silicon_rohf_results_3P_R30.0.jld2",    file_vpol = nothing),
    (element = "Germanium", label = "Germanio (Ge I)", Z = 32,
     file_hf = "germanium_rohf_results_3P_R30.0.jld2",  file_vpol = nothing),
]

# NIST (cm^-1, referidos a 3P_0): [3P0, 3P1, 3P2, 1D2, 1S0]
const NIST = Dict(
    6  => [0.0,  16.4,   43.4, 10192.6, 21648.0],
    14 => [0.0,  77.1,  223.2,  6298.8, 15394.4],
    32 => [0.0, 557.1, 1410.0,  7125.3, 16367.1],
    50 => [0.0, 1691.8, 3427.7, 8613.0, 17162.6],
)

const TERMS = ["^3P_0", "^3P_1", "^3P_2", "^1D_2", "^1S_0"]

function collect_results(cases)
    out = Dict{String,Any}()
    for c in cases
        isfile(c.file_hf) || error("No existe $(c.file_hf). Corre antes el script ROHF de $(c.element).")
        println("### $(c.element): HF+CI desde $(c.file_hf)")
        res_hf = run_toy_ci(c.element, c.file_hf)
        res_hf.angular_ok || error("$(c.element): el test de consistencia angular fallo. " *
                                   "No se genera tabla con numeros que no reproducen 0.24/0.60 F^2.")

        res_vp = res_hf
        if c.file_vpol !== nothing
            isfile(c.file_vpol) || error("No existe $(c.file_vpol).")
            println("### $(c.element): +V_pol desde $(c.file_vpol)")
            res_vp = run_toy_ci(c.element, c.file_vpol)
            res_vp.angular_ok || error("$(c.element) (V_pol): test angular fallido.")
        end
        out[c.element] = (hf = res_hf, vpol = res_vp)
    end
    return out
end

function print_provenance(cases, results)
    println("\n===== PROCEDENCIA (pegar como comentario junto a la tabla) =====")
    for c in cases
        r = results[c.element]
        @printf("%% %-10s HF+CI  : %s | zeta = %.8f Ha | F2 = %.8f Ha | E_corr = %.6e Ha\n",
                c.element, c.file_hf, r.hf.zeta, r.hf.F2, r.hf.E_corr)
        if c.file_vpol !== nothing
            @printf("%% %-10s +Vpol  : %s | zeta = %.8f Ha | F2 = %.8f Ha | E_corr = %.6e Ha\n",
                    c.element, c.file_vpol, r.vpol.zeta, r.vpol.F2, r.vpol.E_corr)
        end
    end
end

function print_latex(cases, results)
    ncol = length(cases)
    println("\n===== LATEX TABLE =====")
    println("\\begin{table}[htbp]")
    println("    \\centering")
    println("    \\resizebox{\\textwidth}{!}{")
    println("    \\begin{tabular}{l|" * join(fill("rrr", ncol), "|") * "}")
    println("        \\toprule")

    header1 = "        "
    for (i, c) in enumerate(cases)
        sep = i == ncol ? "c" : "c|"
        header1 *= "& \\multicolumn{3}{$sep}{\\textbf{$(c.label)}} "
    end
    println(header1 * "\\\\")

    header2 = "        \\textbf{Término} "
    for _ in cases
        header2 *= "& \\textbf{HF+CI} & \\textbf{+V\$_{\\text{pol}}\$} & \\textbf{NIST} "
    end
    println(header2 * "\\\\")
    println("        \\midrule")

    for (i, t) in enumerate(TERMS)
        line = "        \$$t\$ "
        for c in cases
            r = results[c.element]
            line *= @sprintf("& %.1f & %.1f & %.1f ",
                             r.hf.levels[i], r.vpol.levels[i], NIST[c.Z][i])
        end
        println(line * "\\\\")
    end

    println("        \\bottomrule")
    println("    \\end{tabular}")
    println("    }")
    println("    \\caption{Niveles de estructura fina (\$\\text{cm}^{-1}\$) de la secuencia \$np^2\$. " *
            "Se contrasta el modelo de correlación de valencia (HF+CI) con la corrección de " *
            "polarización del core (+\$V_{\\text{pol}}\$) y los valores del NIST.}")
    println("    \\label{tab:niveles_energia}")
    println("\\end{table}")
end

if abspath(PROGRAM_FILE) == @__FILE__
    results = collect_results(CASES)
    print_provenance(CASES, results)
    print_latex(CASES, results)
end
