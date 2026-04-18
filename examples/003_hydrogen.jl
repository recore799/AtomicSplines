using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using LinearAlgebra
using Printf
using Plots

println("\n================================================================")
println("    EXPERIMENT: HYDROGEN p-ORBITALS SUPERPOSITION")
println("================================================================\n")

# ==============================================================================
# A. PHYSICS & SETUP
# ==============================================================================
Z = 1.0
L = 1              # Momento angular l=1 para orbitales p
Target_States = 4  # Visualizaremos 2p, 3p, 4p, 5p

exact_energy(n) = -0.5 * (Z^2) / (n^2)

println("A. Physics Definition:")
println("    System:     Hydrogen (Z = $Z)")
println("    Symmetry:   l = $L (p-orbitals)")
println("    Goal:       Superimpose first $Target_States radial eigenfunctions")

# ==============================================================================
# B. BASIS GENERATION
# ==============================================================================
println("\nB. Discretizing Space...")

R_MAX = 100.0      # Dominio extendido para capturar la cola asintótica del 5p
N_ELEMS = 800      # Mayor resolución para mantener precisión en un dominio grande
ORDER = 7          
GAMMA = 2.5        

basis = generate_basis(R_MAX, N_ELEMS, Val(ORDER), γ=GAMMA)

# ==============================================================================
# C. MATRIX ASSEMBLY & DIAGONALIZATION
# ==============================================================================
println("\nC. Assembling Operators & Diagonalizing...")

ws = init_scf_workspace(basis, Z)
active = 2:(basis.num_splines - 1)

# Construcción del Hamiltoniano Efectivo: Cinética + Coulomb + Centrífuga
centrifugal_factor = L * (L + 1) / 2.0
H = ws.T + ws.V + (centrifugal_factor * ws.V2)

t_solve = @elapsed evals, evecs = eigen(Symmetric(H[active,active]), Symmetric(ws.S[active,active]))
println("    > Diagonalization complete in $(round(t_solve, digits=4)) s")

# ==============================================================================
# D. RESULTS ANALYSIS & VISUALIZATION
# ==============================================================================
println("\nD. Extracting States & Plotting...")

# Definimos una red densa para una evaluación suave de las curvas
r_plot = range(0, 60.0, length=2000) 

# Inicializamos el lienzo principal
p_super = plot(title="Superposición de Orbitales Radiales Hidrogenoides (l=1)", 
               xlabel="Distancia Radial r (a.u.)", 
               ylabel="Amplitud Radial P(r)",
               grid=true, gridalpha=0.3,
               legend=:topright,
               dpi=300)

hline!(p_super, [0], color=:black, lw=1, alpha=0.5, label="")

perm = sortperm(Real.(evals))

# Paleta de colores para distinción clara
colors = [:blue, :red, :green, :purple]

for i in 1:Target_States
    idx = perm[i]
    E_calc = evals[idx]
    
    n = i + L 
    E_ex = exact_energy(n)
    err = abs(E_calc - E_ex)
    
    @printf("  State %dp | E_calc: %12.8f | Error: %.2e\n", n, E_calc, err)
    
    # 1. Extracción y Normalización
    c_full = zeros(Float64, basis.num_splines)
    c_full[active] = evecs[:, idx]
    
    norm_val = sqrt(dot(c_full, ws.S * c_full))
    c_full ./= norm_val
    
    # 2. Evaluación Continua sobre la red r_plot
    psi_vals = evaluate_orbital(basis, c_full, r_plot)
    
    # 3. Convención de Fase (Lóbulo asintótico externo siempre positivo)
    # Buscamos el extremo (máximo o mínimo) que se encuentra más lejos del núcleo.
    # Una heurística robusta es verificar el signo de la función cerca de su valor absoluto máximo global, 
    # asumiendo que el lóbulo externo suele ser el dominante en área, o simplemente invertir si el máximo es menor al mínimo absoluto.
    if maximum(psi_vals) < abs(minimum(psi_vals))
        psi_vals .*= -1.0
    end
    
    # 4. Inyección en la gráfica principal
    label_str = @sprintf("%dp (E=%.4f)", n, E_calc)
    plot!(p_super, r_plot, psi_vals, label=label_str, lw=2.5, color=colors[i], alpha=0.85)
end

display(p_super)
# savefig(p_super, "p_orbitals_superposition.pdf")

println("\n>> Proceso Terminado.")
