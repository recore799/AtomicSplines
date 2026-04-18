using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using LinearAlgebra
using Printf
using Plots

println("\n================================================================")
println("    EXPERIMENT: HYDROGEN ATOM (Spectrum & Orbitals)")
println("================================================================\n")

# ==============================================================================
# A. PHYSICS & SETUP
# ==============================================================================
Z = 1.0
Target_States = 5 

exact_energy(n) = -0.5 * (Z^2) / (n^2)

println("A. Physics Definition:")
println("    System:     Hydrogen (Z = $Z)")
println("    Equation:   H c = E S c")
println("    Goal:       Retrieve first $Target_States radial eigenfunctions")

# ==============================================================================
# B. BASIS GENERATION
# ==============================================================================
println("\nB. Discretizing Space...")

R_MAX = 80.0        
N_ELEMS = 1000        
ORDER = 7          # Spline Order k
GAMMA = 2.5        

# UPDATE 1: Use Val(ORDER) for parametric dispatch
basis = generate_basis(R_MAX, N_ELEMS, Val(ORDER), γ=GAMMA)

println("    > R_max: $R_MAX a.u.")
println("    > Elements: $N_ELEMS (Order $ORDER)")
println("    > Splines: $(basis.num_splines)")

# ==============================================================================
# C. MATRIX ASSEMBLY
# ==============================================================================
println("\nC. Assembling Operators...")

# This builds S, T, V (Banded) and Tensors (unused for hydrogen)
t_build = @elapsed ws = init_scf_workspace(basis, Z)

println("    > Workspace built in $(round(t_build, digits=4)) s")

# 2. Apply Boundary Conditions (Dirichlet)
#    Remove first and last splines to enforce u(0)=0 and u(R)=0
active = 2:(basis.num_splines - 1)

# 3. Construct Hamiltonian
#    ws.T and ws.V are BandedMatrices
H = ws.T + ws.V

# ==============================================================================
# D. SOLVE EIGENPROBLEM
# ==============================================================================
println("\nD. Diagonalizing Hamiltonian...")

t_solve = @elapsed evals, evecs = eigen(Symmetric(H[active,active]), Symmetric(ws.S[active,active]))

println("    > Solved $(length(active)) states in $(round(t_solve, digits=4)) s")

# ==============================================================================
# E. RESULTS ANALYSIS
# ==============================================================================
println("\nE. Spectrum Analysis (First $Target_States states):")
println("----------------------------------------------------------------")
println("State |  Energy (Ha)   |  Exact (Ha)    |  Error")
println("----------------------------------------------------------------")

states_to_plot = []
perm = sortperm(Real.(evals))

for n in 1:Target_States
    idx = perm[n]
    E_calc = evals[idx]
    E_ex = exact_energy(n)
    err = abs(E_calc - E_ex)
    
    @printf("  %ds  | %12.8f   | %12.8f   |  %.2e\n", n, E_calc, E_ex, err)
    
    # Extract Wavefunction Vector
    c_full = zeros(Float64, basis.num_splines)
    c_full[active] = evecs[:, idx]
    
    # Normalize
    norm = sqrt(dot(c_full, ws.S * c_full))
    c_full ./= norm
    
    push!(states_to_plot, (n, E_calc, c_full))
end
println("----------------------------------------------------------------")

# ==============================================================================
# F. VISUALIZATION
# ==============================================================================
println("\nF. Plotting Radial Wavefunctions P(r)...")

r_plot = range(0, 30.0, length=1000) # Zoom in
p = plot(title="Hydrogen Radial Orbitals", 
         xlabel="r (a.u.)", ylabel="P(r)",
         lw=2, legend=:topright)

for (i, (n, E, coeffs)) in enumerate(states_to_plot)
    # Calculate P(r) values
    psi_vals = evaluate_orbital(basis, coeffs, r_plot)
    
    # Phase fix: ensure first lobe is positive for visual consistency
    if maximum(psi_vals) < abs(minimum(psi_vals)); psi_vals .*= -1; end
    
    label_str = @sprintf("%ds (E=%.3f)", n, E)
    # plot!(p, r_plot, psi_vals, label=label_str, color=colors[i], lw=2)
    plot!(p, r_plot, psi_vals, label=label_str, lw=2)
end

hline!(p, [0], color=:black, alpha=0.3, label="")
display(p)
# savefig("hydrogen_spectrum.pdf")

