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
# A. DEFINICIÓN DEL PROBLEMA
#    Resolver la Ec. de Schrödinger Radial para l=0 (simetría esférica)
#    H = -1/2 d²/dr² - Z/r
# ==============================================================================

Z_ATOM = 1.0  # Hidrógeno
exact_energy(n) = -0.5 * (Z_ATOM^2) / (n^2)

println("A. Física del Problema:")
println("   Sistema:   Átomo Hidrogenoide (Z = $Z_ATOM)")
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
# C. ENSAMBLAJE DE HAMILTONIANO
# ==============================================================================
println("\nC. Construyendo Operadores (T, V, S)...")

# AtomicSplines hace el trabajo pesado: integrales cinéticas, nucleares y solapa
T, V, S = assemble_core(basis, Z_ATOM)
H = T + V  # Hamiltoniano total

println("   > Matrices ensambladas. Tamaño: $(size(H))")

# ==============================================================================
# D. SOLUCIÓN (Diagonalización)
# ==============================================================================
println("\nD. Resolviendo Ecuación de Schrödinger...")

# Condiciones de Frontera: ψ(0)=0, ψ(R_max)=0
# Eliminamos el primer y último spline de la base.
inner = 2:(basis.num_splines - 1)

# Convertimos a denso para usar el solver robusto de LAPACK (eigen)
# (Para sistemas grandes usaríamos Arnoldi/Krylov)
H_cut = Matrix(H[inner, inner])
S_cut = Matrix(S[inner, inner])

evals, evecs = eigen(H_cut, S_cut)

println("   > Diagonalización completada.")

# ==============================================================================
# E. ANÁLISIS DE RESULTADOS
# ==============================================================================
println("\nE. Espectroscopía Calculada:")
println("   -------------------------------------------------------------")
println("   Nivel (nl)   Energía Calc. (Ha)    Energía Exacta      Error")
println("   -------------------------------------------------------------")

# Analizamos los primeros 4 estados (1s, 2s, 3s, 4s)
states_to_plot = []

for n in 1:4
    E_calc = evals[n]
    E_exact = exact_energy(n)
    error = abs(E_calc - E_exact)
    
    # Guardamos los vectores para graficar (re-insertando los ceros de frontera)
    psi_vector = vcat(0.0, evecs[:, n], 0.0)
    push!(states_to_plot, (n, E_calc, psi_vector))
    
    @printf("   %ds          % .8f          %.8f          %.1e\n", 
            n, E_calc, E_exact, error)
end
println("   -------------------------------------------------------------")

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
