using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))


using AtomicSplines
using LinearAlgebra
using Printf
using Plots

println("\n================================================================")
println("    EXPERIMENT: B-SPLINE DECOMPOSITION OF HYDROGEN p-ORBITALS")
println("================================================================\n")

# ==============================================================================
# A. PHYSICS & SETUP
# ==============================================================================
Z = 1.0
L = 1              # Momento angular (p-orbital)
Target_States = 3  # Queremos los primeros 3 estados (2p, 3p, 4p)

exact_energy(n) = -0.5 * (Z^2) / (n^2)

# ==============================================================================
# B. BASIS GENERATION
# ==============================================================================
R_MAX = 80.0        
N_ELEMS = 50        
ORDER = 7          
GAMMA = 2.5        

basis = generate_basis(R_MAX, N_ELEMS, Val(ORDER), γ=GAMMA)

# ==============================================================================
# C. MATRIX ASSEMBLY & DIAGONALIZATION
# ==============================================================================
ws = init_scf_workspace(basis, Z)

# Restricción de Condiciones de Frontera de Dirichlet
active = 2:(basis.num_splines - 1)

# Añadimos la barrera centrífuga: l(l+1)/(2r^2)
centrifugal_factor = L * (L + 1) / 2.0
H = ws.T + ws.V + (centrifugal_factor * ws.V2)

evals, evecs = eigen(Symmetric(H[active,active]), Symmetric(ws.S[active,active]))

# ==============================================================================
# D. VISUALIZATION OF BASIS OVERLAP
# ==============================================================================
println("\nGenerating Basis Decomposition Plots for l=$L ...")

r_plot = range(0, 40.0, length=1000) # Ajustado un poco más lejos para los p-orbitales
plots_array = []
perm = sortperm(Real.(evals))

# Iteramos sobre la cantidad de estados solicitados
for i in 1:Target_States
    idx = perm[i]
    E_calc = evals[idx]
    
    # El número cuántico principal real n comienza en L + 1
    n = i + L 
    E_ex = exact_energy(n)
    
    @printf("  State %dp | Calculated: %12.8f Ha | Exact: %12.8f Ha\n", n, E_calc, E_ex)
    
    # Extraer y normalizar el eigenvector
    c_full = zeros(Float64, basis.num_splines)
    c_full[active] = evecs[:, idx]
    norm_val = sqrt(dot(c_full, ws.S * c_full))
    c_full ./= norm_val
    
    # Evaluar la función de onda total P(r)
    psi_total = evaluate_orbital(basis, c_full, r_plot)
    
    # Fijar la fase para consistencia visual (lóbulo exterior dominante positivo)
    if maximum(psi_total) < abs(minimum(psi_total))
        c_full .*= -1
        psi_total .*= -1
    end
    
    title_str = @sprintf("%dp Orbital | E = %.6f Ha", n, E_calc)
    p_state = plot(title=title_str, xlabel=(i == Target_States ? "r (a.u.)" : ""), 
                   ylabel="P(r)", legend=false, grid=true, gridalpha=0.2)
    
    hline!(p_state, [0], color=:black, lw=1, alpha=0.5)
    
    # Dibujar la función de onda total
    plot!(p_state, r_plot, psi_total, color=:black, lw=4, alpha=0.3)
    
    # Dibujar los B-splines individuales ponderados
    temp_coeffs = zeros(Float64, basis.num_splines)
    
    for j in 1:basis.num_splines
        # Filtro de amplitud para optimizar la gráfica
        if abs(c_full[j]) > 1e-3
            fill!(temp_coeffs, 0.0)
            temp_coeffs[j] = c_full[j]
            
            y_spline = evaluate_orbital(basis, temp_coeffs, r_plot)
            plot!(p_state, r_plot, y_spline, lw=1.5, alpha=0.8)
        end
    end
    
    push!(plots_array, p_state)
end

final_plot = plot(plots_array..., layout=(Target_States, 1), size=(800, 300 * Target_States))
display(final_plot)
savefig(final_plot, "hydrogen_p_spline_decomposition.pdf")

