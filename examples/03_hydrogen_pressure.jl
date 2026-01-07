using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using LinearAlgebra
using Plots
using Printf

println("\n================================================================")
println("   EXPERIMENTO: ÁTOMO CONFINADO Y PRESIÓN CUÁNTICA")
println("================================================================\n")

# --- FÍSICA DEL PROBLEMA ---
# Presión Termodinámica: P = - dE / dV
# Volumen Esférico:      V = 4/3 * π * R³
# Regla de la Cadena:    P = - (1 / 4πR²) * (dE / dR)

# Función auxiliar para resolver el H para un radio R dado
function get_energy_and_psi(R_box)
    # 1. Base: Ajustamos los elementos para mantener densidad constante
    # Si R es pequeño, necesitamos menos elementos; si es grande, más.
    # Densidad aprox: 2 elementos por u.a.
    n_elems = max(20, Int(round(2.0 * R_box)))
    order = 5
    
    basis = generate_basis(R_box, n_elems, order, γ=1.5) # Gamma bajo para cajas pequeñas
    
    # 2. Ensamblaje (Z=1 Hidrógeno)
    T, V, S = assemble_core(basis, 1.0)
    
    # 3. Solver
    H = T + V
    
    # Condiciones de frontera (Caja dura en R_box)
    inner = 2:(basis.num_splines - 1)
    
    vals, vecs = eigen(Matrix(H[inner, inner]), Matrix(S[inner, inner]))
    
    # Retornamos Energía base y coeficientes de la función de onda
    return vals[1], vcat(0.0, vecs[:, 1], 0.0), basis
end

# ==============================================================================
# 1. EXPERIMENTO VISUAL: Función de Onda vs Confinamiento
# ==============================================================================
println("1. Calculando funciones de onda para diferentes radios...")

radii_visual = [1.1, 1.835, 2.5, 10.0] # Radios de confinamiento (u.a.)
densities = []

p1 = plot(title="Densidad Radial de Probabilidad |P(r)|²",
          xlabel="Distancia r (u.a.)", ylabel="Densidad",
          xlims=(0, 6), legend=:topright)

for R in radii_visual
    E, coeffs, basis = get_energy_and_psi(R)
    
    # Graficamos en una rejilla común para comparar
    r_grid = range(0, R, length=300)
    psi_vals = [eval_expansion(coeffs, basis, r) for r in r_grid]
    
    # Normalización Numérica simple (suma de Riemann) para que el área sea ~1
    dr = r_grid[2] - r_grid[1]
    norm_factor = sum(psi_vals.^2) * dr
    prob_density = (psi_vals.^2) ./ norm_factor
    
    @printf("   > Radio de Caja: %4.1f u.a. | Energía: %7.4f Ha\n", R, E)
    
    plot!(p1, r_grid, prob_density, lw=2, 
          label="R_caja = $(R) (E=$(@sprintf("%.3f", E)))", fill=(0, 0.2))
end

# Dibujar una pared vertical simulada para R=2 para visualizar el choque
vline!(p1, [2.0], color=:black, ls=:dash, label="Pared R=2")


# ==============================================================================
# 2. CÁLCULO DE PRESIÓN (Curva P vs R)
# ==============================================================================
println("\n2. Calculando curva de Presión Cuántica...")

R_scan = range(1.5, 8.0, length=30) # Barrido de radios
Energies = zeros(length(R_scan))
Pressures = zeros(length(R_scan))

# Diferencia finita para la derivada dE/dR
dR = 1e-4 

for (i, R) in enumerate(R_scan)
    E_current, _, _ = get_energy_and_psi(R)
    E_plus, _, _    = get_energy_and_psi(R + dR)
    
    # Derivada Numérica
    dEdR = (E_plus - E_current) / dR
    
    # Fórmula de Presión
    # P = - (1 / Area) * (dE / dR)
    area = 4 * π * R^2
    P = - dEdR / area
    
    Energies[i] = E_current
    Pressures[i] = P
end

# Gráfica de Energía
p2 = plot(R_scan, Energies, lw=2, color=:blue, legend=false,
          xlabel="Radio de Confinamiento (u.a.)", ylabel="Energía (Ha)", title="Energía vs Radio")
hline!(p2, [-0.5], color=:red, ls=:dash, label="Límite Libre (-0.5)")

# Gráfica de Presión
p3 = plot(R_scan, Pressures, lw=2, color=:red, legend=false,
          xlabel="Radio de Confinamiento (u.a.)", ylabel="Presión (u.a.)",
          title="Presión Cuántica")

# Combinar todo
l = @layout [a; [b c]]
final_plot = plot(p1, p2, p3, layout=l, size=(800, 800))
display(final_plot)

