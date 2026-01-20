# 1. Configuración del Entorno
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

# 2. Carga de Herramientas
using AtomicSplines
using LinearAlgebra
using FastGaussQuadrature
using Plots
using Printf

println("\n================================================================")
println("   EXPERIMENTO: SOLVER DE POISSON 1D (FEM B-SPLINES)")
println("================================================================\n")

# ==============================================================================
# A. DEFINICIÓN DEL PROBLEMA
# ==============================================================================
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
N_ELEMS = 20
ORDER = 4 # B-splines cúbicos

basis = generate_basis(R_MAX, N_ELEMS, ORDER, γ=1.0)

println("   > Parámetros: R_max=$R_MAX, Elementos=$N_ELEMS, Orden=$ORDER")
println("   > Total Splines: $(basis.num_splines)")
println("   > Vector de Nudos (Knots):")
# Mostramos los nudos formateados para ver la multiplicidad en los bordes
pretty_knots = [(@sprintf("%.2f", k)) for k in basis.knots]
println("     [$(join(pretty_knots, ", "))]")
println("     (Nota: La repetición en 0.00 y 1.00 se llama 'clamping')\n")

# ==============================================================================
# C. ENSAMBLAJE DE MATRICES
# ==============================================================================
println("C. Ensamblando Sistema Matricial...")

n = basis.num_splines
K = zeros(n, n)
F = zeros(n)
gl_p, gl_w = gausslegendre(basis.order + 2)

# Contador para verificar dispersión
non_zeros = 0

for i in 1:n
    for j in 1:n
        if abs(i - j) >= basis.order
            continue 
        end
        
        val_k = 0.0
        start_knot = max(i, j)
        end_knot   = min(i + basis.order, j + basis.order)
        
        for e in start_knot:(end_knot - 1)
            a = basis.knots[e]; b = basis.knots[e+1]
            if a == b; continue; end
            mid = (a + b) / 2; scale = (b - a) / 2
            
            for q in 1:length(gl_p)
                x = mid + scale * gl_p[q]
                w = scale * gl_w[q]
                val_k += w * d_bspline(i, basis.order, x, basis.knots) * d_bspline(j, basis.order, x, basis.knots)
            end
        end
        K[i,j] = val_k
        if abs(val_k) > 1e-12; global non_zeros += 1; end
    end
    
    val_f = 0.0
    start_knot = i; end_knot = i + basis.order
    for e in start_knot:(end_knot - 1)
        a = basis.knots[e]; b = basis.knots[e+1]
        if a == b; continue; end
        mid = (a + b) / 2; scale = (b - a) / 2
        for q in 1:length(gl_p)
            x = mid + scale * gl_p[q]; w = scale * gl_w[q]
            val_f += w * bspline(i, basis.order, x, basis.knots) * source_f(x)
        end
    end
    F[i] = val_f
end

sparsity = (non_zeros / (n^2)) * 100
println("   > Matriz K ensamblada ($n x $n)")
println("   > Elementos no nulos: $non_zeros")
@printf("   > Dispersión (Sparsity): %.1f%% llena (¡El resto son ceros!)\n\n", sparsity)

# ==============================================================================
# D. SOLUCIÓN
# ==============================================================================
println("D. Resolviendo K * c = F...")

inner_indices = 2:(n-1)
coeffs = zeros(n)
coeffs[inner_indices] = K[inner_indices, inner_indices] \ F[inner_indices]

println("   > Sistema lineal resuelto.\n")

# ==============================================================================
# 4. ANÁLISIS DE RESULTADOS
# ==============================================================================
println("E. Análisis de Resultados...")

x_plot = range(0, R_MAX, length=200)
u_exact = exact_u.(x_plot)
u_fem   = zeros(length(x_plot))
weighted_splines = zeros(length(x_plot), n)

for (k, x) in enumerate(x_plot)
    val_sum = 0.0
    for i in 1:n
        contribution = coeffs[i] * bspline(i, basis.order, x, basis.knots)
        weighted_splines[k, i] = contribution
        val_sum += contribution
    end
    u_fem[k] = val_sum
end

# Cálculo de Error L2
error_l2 = norm(u_fem - u_exact) / norm(u_exact)
@printf("   > Error Relativo L2: %.2e  <-- (Debería ser muy bajo)\n", error_l2)

# Gráficas
p1 = plot(x_plot, weighted_splines, legend=false, color=:blue, alpha=0.5, 
          title="Contribución de cada B-Spline")
plot!(p1, x_plot, u_fem, color=:red, lw=3, label="Suma FEM")

p2 = plot(x_plot, u_exact, label="Exacta", color=:black, ls=:dash, lw=2, 
          title="Comparación Final")
plot!(p2, x_plot, u_fem, label="B-Spline FEM", color=:red, alpha=0.7, lw=2)

final_plot = plot(p1, p2, layout=(2,1), size=(600, 800))
display(final_plot)

println("\n>> ¡Experimento completado con éxito!")
