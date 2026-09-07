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
