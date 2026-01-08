using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using LinearAlgebra
using Plots
using Printf

println("\n================================================================")
println("   EXPERIMENTO: ÁTOMO CONFINADO Y PRESIÓN (API V2)")
println("================================================================\n")

# --- LÓGICA DE FÍSICA ---
function get_confinement_data(R_box)
    # 1. Definimos el átomo (Hidrógeno 1s)
    h_atom = Atom(1.0, [Orbital(1, 0, 1.0)])
    
    # 2. Base dinámica según el radio
    n_elems = max(25, Int(round(2.5 * R_box)))
    basis = generate_basis(R_box, n_elems, 5, γ=1.2)
    
    # 3. Uso de la nueva API
    # solve_orbital! se encarga de: assemble_hamiltonian -> solve_eigen -> update orbital
    solve_orbital!(h_atom.orbitals[1], h_atom, basis)
    
    return h_atom.orbitals[1], basis
end

# ==============================================================================
# 1. VISUALIZACIÓN DE DENSIDADES
# ==============================================================================
println("1. Generando densidades radiales para diferentes confinamientos...")

radii_visual = [1.1, 1.835, 2.5, 5.0]
p1 = plot(title="Densidad Radial |P(r)|² bajo Confinamiento",
          xlabel="Distancia r (u.a.)", ylabel="Densidad",
          xlims=(0, 6), legend=:topright)

for R in radii_visual
    orb, basis = get_confinement_data(R)
    
    r_grid = range(0, R, length=400)
    # eval_expansion usa orb.coeffs directamente
    rho = [eval_expansion(orb.coeffs, basis, r)^2 for r in r_grid]
    
    # Normalización para asegurar unidad en la integral
    dr = r_grid[2] - r_grid[1]
    rho ./= (sum(rho) * dr)
    
    plot!(p1, r_grid, rho, lw=2, fill=(0, 0.1),
          label="R=$R (E=$(@sprintf("%.3f", orb.energy)))")
end

# ==============================================================================
# 2. CÁLCULO DE PRESIÓN TERMODINÁMICA
# ==============================================================================
println("\n2. Calculando curva de Presión Cuántica P = -dE/dV...")

R_scan = range(1.2, 7.0, length=35)
Energies = zeros(length(R_scan))
Pressures = zeros(length(R_scan))
dR = 1e-4

for (i, R) in enumerate(R_scan)
    orb_now, _ = get_confinement_data(R)
    orb_plus, _ = get_confinement_data(R + dR)
    
    # dE/dR numérico
    dEdR = (orb_plus.energy - orb_now.energy) / dR
    
    # P = -(1/4πR²) * (dE/dR)
    P = -dEdR / (4 * π * R^2)
    
    Energies[i] = orb_now.energy
    Pressures[i] = P
end

# --- GRÁFICAS ---
p2 = plot(R_scan, Energies, lw=2, color=:blue, legend=false,
          xlabel="Radio R (u.a.)", ylabel="Energía (Ha)", title="Energía vs Confinamiento")
hline!(p2, [-0.5], color=:black, ls=:dash)

p3 = plot(R_scan, Pressures, lw=2, color=:red, legend=false,
          xlabel="Radio R (u.a.)", ylabel="Presión (u.a.)", title="Presión sobre la Pared")

l = @layout [a; [b c]]
final_plot = plot(p1, p2, p3, layout=l, size=(900, 850))
display(final_plot)

println("\n>> ¡Éxito! El rediseño procesó $(length(R_scan)*2 + length(radii_visual)) diagonalizaciones correctamente.")
