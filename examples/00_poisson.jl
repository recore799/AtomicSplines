# 1. Configuración del Entorno
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using LinearAlgebra
using FastGaussQuadrature
using Plots
using Printf

println("\n================================================================")
println("   EXPERIMENTO: ECUACIÓN DE POISSON 1D (FEM B-SPLINES / API V2)")
println("================================================================\n")

# ==============================================================================
# A. DEFINICIÓN DEL PROBLEMA
# ==============================================================================
# Problema: -u''(x) = f(x)
exact_u(x) = sin(2*π*x)
source_f(x) = (4*π^2) * sin(2*π*x)

println("A. Problema Definido:")
println("   Ecuación: -u''(x) = f(x)")
println("   Dominio:  [0, 1]")
println("   BCs:      u(0) = u(1) = 0\n")

# ==============================================================================
# B. GENERACIÓN DE LA BASE
# ==============================================================================
println("B. Generando base B-Spline...")

R_MAX = 1.0
N_ELEMS = 5
ORDER = 4     # Cúbicos
basis = generate_basis(R_MAX, N_ELEMS, ORDER, γ=1.0) # γ=1.0 es malla uniforme

println("   > Base: $N_ELEMS elementos, Orden $ORDER")
println("   > Total Splines: $(basis.num_splines)")

# ==============================================================================
# C. ENSAMBLAJE DE MATRICES
# ==============================================================================
println("C. Ensamblando Sistema Matricial...")

# 1. Matriz de Rigidez (Stiffness Matrix) K
#    Sabemos que T_ij = 0.5 * integral(B_i' * B_j')
#    La ecuación de Poisson requiere K_ij = integral(B_i' * B_j')
#    Por lo tanto: K = 2 * T
print("   > Reutilizando matriz cinética para K... ")
T, _, S = assemble_core(basis, 0.0) # Z=0 porque no hay núcleo aquí
K = 2.0 .* T
println("Listo.")

# 2. Vector de Carga (Load Vector) F
#    F_i = integral(B_i(x) * f(x) dx)
#    Esto requiere un ensamblaje manual porque f(x) es arbitraria
println("   > Ensamblando vector de carga F (Integración manual)...")

n = basis.num_splines
F = zeros(n)
k = basis.order

# Buffers pre-asignados (API V2)
vals = zeros(Float64, k)
derivs = zeros(Float64, k) # No se usa aquí pero el kernel lo pide
gl_p, gl_w = gausslegendre(k + 2)

# Iteramos sobre ELEMENTOS (intervalos de nudos)
for i in 1:(length(basis.knots)-1)
    t_a = basis.knots[i]; t_b = basis.knots[i+1]
    if t_a == t_b; continue; end # Intervalo vacío (multiplicidad)

    mid = (t_a + t_b)/2
    scale = (t_b - t_a)/2
    first_global = i - k + 1

    for q in 1:length(gl_p)
        x = scale * gl_p[q] + mid
        w = scale * gl_w[q]
        fx = source_f(x)

        # Llamada al Kernel de bajo nivel (solo valores)
        # Nota: Usamos la función interna eval_bspline_kernel! que ya conoces
        AtomicSplines.eval_bspline_kernel!(vals, derivs, Val(true), Val(false), 
                                           i, k, x, basis.knots)

        # "Scatter" a los índices globales
        for local_idx in 1:k
            global_idx = first_global + local_idx - 1
            if global_idx >= 1 && global_idx <= n
                F[global_idx] += w * vals[local_idx] * fx
            end
        end
    end
end
println("   > Vector F completado.")

# ==============================================================================
# D. SOLUCIÓN
# ==============================================================================
println("\nD. Resolviendo K * c = F...")

# Condiciones de frontera: u(0)=0, u(1)=0
# En B-Splines (con nudos repetidos en bordes), esto implica c_1 = 0 y c_n = 0
active = 2:(n-1)

coeffs = zeros(n)
# Resolvemos el sub-sistema
coeffs[active] = K[active, active] \ F[active]

println("   > Sistema lineal resuelto.")

# ==============================================================================
# E. ANÁLISIS DE RESULTADOS
# ==============================================================================
println("E. Análisis de Resultados...")

# 1. Definir red de puntos para graficar
x_plot = collect(range(0, R_MAX, length=200))

# 2. Solución Exacta
u_ex = exact_u.(x_plot)

# 3. Solución FEM (¡Usando tu nueva función vectorizada!)
#    Esto reemplaza el bucle doble manual. Es más rápido y limpio.
u_fem = evaluate_orbital(basis, coeffs, x_plot)

# 4. Cálculo de Error L2
#    Norma Euclidiana del vector diferencia
err = norm(u_fem - u_ex) / norm(u_ex)
@printf("   > Error Relativo L2: %.2e\n", err)

# --- GRÁFICAS ---

# Gráfica 1: Comparación Directa
p1 = plot(x_plot, u_ex, label="Exacta - sin(2πx)", color=:black, ls=:dash, lw=2,
          title="Solución de Poisson 1D", xlabel="x", ylabel="u(x)")
plot!(p1, x_plot, u_fem, label="B-Spline FEM (k=$ORDER)", color=:red, alpha=0.7, lw=3)

# Gráfica 2: Descomposición de la Base
# Queremos ver cómo cada spline individual contribuye a la suma final.
p2 = plot(title="Contribución por Spline", legend=false, xlabel="x")

# Dibujamos la suma total como referencia (sombra gris)
plot!(p2, x_plot, u_fem, color=:black, lw=0, fill=(0, 0.1, :black))

# Dibujamos cada spline individual ponderado por su coeficiente c_i
# Truco: Reutilizamos evaluate_orbital pasando un vector con un solo '1'
temp_coeffs = zeros(basis.num_splines)

for i in 1:basis.num_splines
    # Solo graficamos si el coeficiente es significativo (para no llenar de líneas planas)
    if abs(coeffs[i]) > 1e-3
        fill!(temp_coeffs, 0.0)
        temp_coeffs[i] = coeffs[i] # Solo activamos el spline 'i'
        
        # Evaluamos solo este spline
        y_spline = evaluate_orbital(basis, temp_coeffs, x_plot)
        
        plot!(p2, x_plot, y_spline, lw=1.5, alpha=0.8)
    end
end

final_plot = plot(p1, p2, layout=(2,1), size=(600, 800))
display(final_plot)
savefig(final_plot, "poisson_1d_solution.png")

println("\n>> Final")
