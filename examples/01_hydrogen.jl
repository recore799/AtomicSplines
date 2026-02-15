using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using LinearAlgebra
using Printf
using Plots
using BandedMatrices

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
N_ELEMS = 100        
ORDER = 6          # Spline Order k
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

# Convert BandedMatrix to Dense Matrix for standard 'eigen'
H_active = Matrix(H[active, active])
S_active = Matrix(ws.S[active, active])

t_solve = @elapsed evals, evecs = eigen(H_active, S_active)

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


# println("\n================================================================")
# println("   EXPERIMENT: HYDROGEN ATOM (Spectrum & Orbitals)")
# println("================================================================\n")

# # ==============================================================================
# # A. PHYSICS & SETUP
# # ==============================================================================
# Z = 1.0
# Target_States = 5 # We want 1s, 2s, 3s, 4s

# exact_energy(n) = -0.5 * (Z^2) / (n^2)

# println("A. Physics Definition:")
# println("   System:    Hydrogen (Z = $Z)")
# println("   Equation:  H c = E S c")
# println("   Goal:      Retrieve first $Target_States radial eigenfunctions")

# # ==============================================================================
# # B. BASIS GENERATION
# # ==============================================================================
# println("\nB. Discretizing Space...")

# # Overkill basis for ground state, but required for excited states
# R_MAX = 100.0       
# N_ELEMS = 100       
# ORDER = 6          
# GAMMA = 2.5        

# basis = generate_basis(R_MAX, N_ELEMS, ORDER, γ=GAMMA)

# println("   > R_max: $R_MAX a.u.")
# println("   > Elements: $N_ELEMS (Order $ORDER)")
# println("   > Splines: $(basis.num_splines)")

# # ==============================================================================
# # C. MATRIX ASSEMBLY
# # ==============================================================================
# println("\nC. Assembling Operators...")

# # 1. Build Matrices
# #    For Hydrogen, we only need the Core Hamiltonian. 
# #    There is no SCF loop because there is no electron-electron repulsion.
# t_build = @elapsed T, V_nuc, S = assemble_core(basis, Z)
# println("   > Matrices built in $(round(t_build, digits=4)) s")

# # 2. Apply Boundary Conditions (Dirichlet)
# #    Remove first and last splines to enforce u(0)=0 and u(R)=0
# active = 2:(basis.num_splines - 1)

# # 3. Construct Hamiltonian
# H = T + V_nuc

# # ==============================================================================
# # D. SOLVE EIGENPROBLEM
# # ==============================================================================
# println("\nD. Diagonalizing Hamiltonian...")

# t_solve = @elapsed evals, evecs = eigen(H[active, active], S[active, active])

# println("   > Solved $(length(active)) states in $(round(t_solve, digits=4)) s")

# # ==============================================================================
# # E. RESULTS ANALYSIS
# # ==============================================================================
# println("\nE. Spectrum Analysis (First $Target_States states):")
# println("----------------------------------------------------------------")
# println("State |  Energy (Ha)   |  Exact (Ha)    |  Error")
# println("----------------------------------------------------------------")

# states_to_plot = []

# # Indices in the sorted eigen-spectrum
# perm = sortperm(Real.(evals))

# for n in 1:Target_States
#     # Get the n-th root
#     idx = perm[n]
#     E_calc = evals[idx]
    
#     # Compare to exact energy
#     E_ex = exact_energy(n)
#     err = abs(E_calc - E_ex)
    
#     @printf("  %ds  | %12.8f   | %12.8f   |  %.2e\n", n, E_calc, E_ex, err)
    
#     # Extract Wavefunction Vector
#     c_full = zeros(Float64, basis.num_splines)
#     c_full[active] = evecs[:, idx]
    
#     # Normalize (just in case)
#     norm = sqrt(dot(c_full, S * c_full))
#     c_full ./= norm
    
#     push!(states_to_plot, (n, E_calc, c_full))
# end
# println("----------------------------------------------------------------")

# # ==============================================================================
# # F. VISUALIZATION
# # ==============================================================================
# println("\nF. Plotting Radial Wavefunctions P(r)...")

# r_plot = range(0, 30.0, length=1000) # Zoom in
# p = plot(title="Hydrogen Radial Orbitals", 
#          xlabel="r (a.u.)", ylabel="P(r)",
#          lw=2, legend=:topright)

# for (i, (n, E, coeffs)) in enumerate(states_to_plot)
#     # Calculate P(r) values
#     psi_vals = evaluate_orbital(basis, coeffs, r_plot)
    
#     # Phase fix: ensure first lobe is positive for visual consistency
#     if maximum(psi_vals) < abs(minimum(psi_vals)); psi_vals .*= -1; end
    
#     label_str = @sprintf("%ds (E=%.3f)", n, E)
#     # plot!(p, r_plot, psi_vals, label=label_str, color=colors[i], lw=2)
#     plot!(p, r_plot, psi_vals, label=label_str, lw=2)
# end

# hline!(p, [0], color=:black, alpha=0.3, label="")
# display(p)
# savefig("hydrogen_spectrum.pdf")
