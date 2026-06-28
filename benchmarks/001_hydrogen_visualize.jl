# 1. Configuración del Entorno
using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using AtomicSplines
using LinearAlgebra
using Printf
using Plots

println("\n================================================================")
println("    EXPERIMENT: B-SPLINE DECOMPOSITION OF HYDROGEN ORBITALS")
println("================================================================\n")

# ==============================================================================
# A. PHYSICS & SETUP
# ==============================================================================
Z = 1.0
Target_States = 3 

exact_energy(n) = -0.5 * (Z^2) / (n^2)

# ==============================================================================
# B. BASIS GENERATION
# ==============================================================================
R_MAX = 60.0        
N_ELEMS = 20        
ORDER = 7          
GAMMA = 2.5        

basis = generate_basis(R_MAX, N_ELEMS, Val(ORDER), γ=GAMMA)

# ==============================================================================
# C. MATRIX ASSEMBLY & DIAGONALIZATION
# ==============================================================================
ws = init_scf_workspace(basis, Z)

active = 2:(basis.num_splines - 1)
H = ws.T + ws.V

evals, evecs = eigen(Symmetric(H[active,active]), Symmetric(ws.S[active,active]))

# ==============================================================================
# D. VISUALIZATION OF BASIS OVERLAP
# ==============================================================================
println("\nGenerating Basis Decomposition Plots...")

r_plot = range(0, 30.0, length=1000) 
plots_array = []
perm = sortperm(Real.(evals))

for n in 1:Target_States
    idx = perm[n]
    E_calc = evals[idx]
    
    # Extract and normalize the eigenvector
    c_full = zeros(Float64, basis.num_splines)
    c_full[active] = evecs[:, idx]
    norm_val = sqrt(dot(c_full, ws.S * c_full))
    c_full ./= norm_val
    
    # Evaluate total wavefunction
    psi_total = evaluate_orbital(basis, c_full, r_plot)
    
    # Phase fix: ensure the dominant outer lobe is positive for standard convention
    if maximum(psi_total) < abs(minimum(psi_total))
        c_full .*= -1
        psi_total .*= -1
    end
    
    # Initialize the subplot for this specific state
    title_str = @sprintf("%ds Orbital | E = %.6f Ha", n, E_calc)
    p_state = plot(title=title_str, xlabel=(n == Target_States ? "r (a.u.)" : ""), 
                   ylabel="P(r)", legend=false, grid=true, gridalpha=0.2)
    
    # Draw the zero-line reference
    hline!(p_state, [0], color=:black, lw=1, alpha=0.5)
    
    # Plot the total wavefunction as a thick, semi-transparent backdrop
    plot!(p_state, r_plot, psi_total, color=:black, lw=4, alpha=0.3)
    
    # Plot individual weighted B-splines
    temp_coeffs = zeros(Float64, basis.num_splines)
    
    for i in 1:basis.num_splines
        # Filter out negligible splines to optimize plotting speed and visual clarity
        if abs(c_full[i]) > 1e-3
            fill!(temp_coeffs, 0.0)
            temp_coeffs[i] = c_full[i]
            
            y_spline = evaluate_orbital(basis, temp_coeffs, r_plot)
            
            # Use alternating colors or a single color map to distinguish splines
            plot!(p_state, r_plot, y_spline, lw=1.5, alpha=0.8)
        end
    end
    
    push!(plots_array, p_state)
end

# Combine all subplots into a single vertical layout
final_plot = plot(plots_array..., layout=(Target_States, 1), size=(800, 300 * Target_States))
display(final_plot)
savefig(final_plot, "hydrogen_spline_decomposition.pdf")

println("Done. Graphs displayed.")
