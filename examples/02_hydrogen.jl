using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using LinearAlgebra
using Printf
using Plots

println("\n================================================================")
println("   EXPERIMENTO: ÁTOMO DE HIDRÓGENO (Espectro y Orbitales)")
println("================================================================\n")

# ==============================================================================
# A. DEFINICIÓN DEL ÁTOMO
# ==============================================================================

my_orbitals = [
    Orbital(1,0,1.0)
    Orbital(2,0,0.0)
    Orbital(3,0,0.0)
    Orbital(4,0,0.0)
]

hydrogen = Atom(1.0, my_orbitals)

exact_energy(n) = -0.5 * (hydrogen.Z^2) / (n^2)

println("A. Física del Problema:")
println("   Sistema:   Átomo Hidrogenoide (Z = $(hydrogen.Z))")
println("   Ecuación:  H ψ = E S ψ  (Problema de valores propios generalizado)")
println("   Teoría:    E_n = -0.5 / n^2 Ha")

# ==============================================================================
# B. GENERACIÓN DE LA BASE (Malla Exponencial)
# ==============================================================================
println("\nB. Discretización del Espacio...")

# Parámetros ajustados para capturar tanto el pico en r=0 como la cola exponencial
R_MAX = 60.0       # Caja grande para contener estados excitados difusos
N_ELEMS = 60       # Pocos elementos, pero bien distribuidos
ORDER = 6          # Alto orden para precisión espectral
GAMMA = 2.5        # Factor de concentración cerca del núcleo

basis = generate_basis(R_MAX, N_ELEMS, ORDER, γ=GAMMA)

println("   > R_max: $R_MAX u.a.")
println("   > Elementos: $N_ELEMS (Orden $ORDER)")
println("   > Distribución: Exponencial (γ=$GAMMA)")
println("   > Total Splines: $(basis.num_splines)")

# ==============================================================================
# D. SOLUCIÓN
# ==============================================================================
println("\nD. Resolviendo cada orbital definido en el Átomo...")

states_to_plot = []

for orb in hydrogen.orbitals
    # Esta funcion hace TODO: ensambla H (con su l), resuelve y guarda E y coeffs
    solve_orbital!(orb, hydrogen, basis)

    # Calculamos error para el print
    E_ex = exact_energy(orb.n)
    err = abs(orb.energy - E_ex)

    @printf("    Orbital %d%s: E_calc = %.8f | E_exact = %.8f | Err = %.1e\n",
            orb.n, orb.l == 0 ? "s" : "p", orb.energy, E_ex, err)

    # Guardamos para la gráfica
    push!(states_to_plot, (orb.n, orb.energy, orb.coeffs))
end


println("   > Diagonalización completada.")

# ==============================================================================
# F. VISUALIZACIÓN DE ORBITALES
# ==============================================================================
println("\nF. Graficando Funciones de Onda Radiales P(r)...")

# Graficamos solo los primeros 20 a.u. donde ocurre la física interesante
r_plot = range(0, 20.0, length=500)
p = plot(title="Orbitales Radiales del Hidrógeno (FEM B-Spline)", 
         xlabel="Distancia r (u.a.)", ylabel="P(r) = r R(r)",
         lw=2, legend=:bottomright)

# Añadimos cada estado al gráfico
colors = [:blue, :red, :green, :purple]

for (i, (n, E, coeffs)) in enumerate(states_to_plot)
    # Reconstruimos la función continua a partir de los coeficientes B-spline
    psi_vals = [eval_expansion(coeffs, basis, r) for r in r_plot]
    
    # Normalización simple para visualización (pico máximo = ±1)
    psi_vals /= maximum(abs.(psi_vals))
    # Ajuste de fase (para que empiecen positivos)
    if psi_vals[2] < 0; psi_vals .*= -1; end
    
    label_str = @sprintf("%ds (E=%.3f)", n, E)
    plot!(p, r_plot, psi_vals, label=label_str, color=colors[i], lw=2)
end

# Línea base cero
hline!(p, [0], color=:black, alpha=0.3, label="")

display(p)
