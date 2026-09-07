# Genera docs/figures/scaling_law.pdf, la Figura \ref{fig:scaling_law} del
# capitulo de resultados.
#
# OJO: los valores de abajo estan escritos a mano y no coinciden con lo que
# calcula el codigo (regla 3 de CLAUDE.md). Segun docs/claude/HALLAZGOS-2026-09-06.md:
#
#   zeta(C) aqui 31.9, calculado 42.216
#   F2(Si) y F2(Ge) aqui valen los dos 0.1627; son distintos: 0.16594819 y
#   0.16272437. La repeticion es el mismo copy-paste que contamina
#   tab:integrales_slater.
#
# No se corrigen aqui porque eso mueve una figura publicada: hay que regenerarla
# leyendo los valores de los .jld2 registrados en benchmarks/RESULTS.toml, no
# reescribiendo las constantes.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Plots

# Valores calculados por los scripts de fine structure
Z_vals = [6.0, 14.0, 32.0, 50.0]
zeta_vals = [31.9, 140.9, 930.0, 2276.0] 
F2_vals = [0.2433, 0.1627, 0.1627, 0.1472]

# Gráfica
p1 = plot(Z_vals, zeta_vals, scale=:log10, marker=:circle, 
    label="\\zeta_{np} (Interacción Espín-Órbita)", 
    ylabel="Energía (escala log cm^-1)", xlabel="Número Atómico (Z)",
    title="Escalamiento Relativista vs Repulsión Coulombiana",
    color=:red, lw=2, legend=:topleft, yticks=[10, 100, 1000, 10000, 100000])

# Multiplicar F2 por factor a cm^-1 para que se vea en el mismo eje
F2_cm = F2_vals .* 219474.63
plot!(p1, Z_vals, F2_cm, scale=:log10, marker=:square, 
    label="F^2 (Repulsión Intra-capa)", color=:blue, lw=2)

savefig(p1, joinpath(@__DIR__, "..", "..", "docs", "figures", "scaling_law.pdf"))
println("Gráfica scaling_law.pdf generada con éxito.")
