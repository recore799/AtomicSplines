using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using LinearAlgebra
using Plots
using Printf

println("\n================================================================")
println("   EXPERIMENT: CONFINED ATOM & QUANTUM PRESSURE (New API)")
println("================================================================\n")

"""
Helper function to solve the Ground State of Hydrogen in a box of radius R.
Returns: Energy, Coefficients, Basis
"""
function solve_confined_hydrogen(R_box)
    # 1. Define Basis scaled to the box size
    #    We adjust N_elements to maintain resolution as the box grows/shrinks
    n_elems = max(25, Int(round(2.5 * R_box)))
    ORDER = 5
    GAMMA = 1.2 # Less aggressive clustering for confined boxes (we need resolution at the wall)
    
    basis = generate_basis(R_box, n_elems, ORDER, γ=GAMMA)
    
    # 2. Assemble Operators (Z=1.0)
    T, V, S = assemble_core(basis, 1.0)
    
    # 3. Apply Hard-Wall Boundary Conditions (Dirichlet)
    #    c[1]=0 -> u(0)=0 (Regularity at origin)
    #    c[N]=0 -> u(R_box)=0 (Infinite Potential Wall)
    active = 2:(basis.num_splines - 1)
    
    H = T + V
    
    # 4. Solve Eigenproblem
    #    We only need the ground state (lowest eigenvalue)
    evals, evecs = eigen(H[active, active], S[active, active])
    
    # Sort to ensure we get the ground state
    perm = sortperm(Real.(evals))
    E_ground = evals[perm[1]]
    
    # Reconstruct full coefficient vector
    coeffs = zeros(Float64, basis.num_splines)
    coeffs[active] = evecs[:, perm[1]]
    
    # Normalize (L2 norm)
    norm = sqrt(dot(coeffs, S * coeffs))
    coeffs ./= norm
    
    return E_ground, coeffs, basis
end

# ==============================================================================
# 1. DENSITY VISUALIZATION
# ==============================================================================
println("1. Generating radial densities for different confinement radii...")

radii_visual = [1.1, 1.835, 2.5, 5.0]

p1 = plot(title="Radial Density |P(r)|² under Confinement",
          xlabel="Distance r (a.u.)", ylabel="Density",
          xlims=(0, 6), legend=:topright)

for R in radii_visual
    E, coeffs, basis = solve_confined_hydrogen(R)
    
    # Generate grid up to the wall
    r_grid = collect(range(0, R, length=400))
    
    # NEW API: Fast vectorized evaluation
    psi = evaluate_orbital(basis, coeffs, r_grid)
    rho = psi .^ 2
    
    # Normalize density for plotting (Area = 1)
    dr = r_grid[2] - r_grid[1]
    rho ./= (sum(rho) * dr)
    
    plot!(p1, r_grid, rho, lw=2, fill=(0, 0.1),
          label="R=$R (E=$(@sprintf("%.3f", E)))")
end

# ==============================================================================
# 2. QUANTUM PRESSURE CALCULATION (P = -dE/dV)
# ==============================================================================
println("\n2. Calculating Quantum Pressure Curve...")

R_scan = range(1.2, 7.0, length=35)
Energies = zeros(length(R_scan))
Pressures = zeros(length(R_scan))

# Small delta for numerical derivative
dR = 1e-4

for (i, R) in enumerate(R_scan)
    # Solve at R
    E_now, _, _ = solve_confined_hydrogen(R)
    
    # Solve at R + dR
    E_plus, _, _ = solve_confined_hydrogen(R + dR)
    
    # Numerical Derivative: dE/dR
    dEdR = (E_plus - E_now) / dR
    
    # Thermodynamics: Pressure = -dE/dV
    # Volume V = (4/3) * pi * R^3
    # dV/dR = 4 * pi * R^2
    # P = -(dE/dR) / (dV/dR)
    P = -dEdR / (4 * π * R^2)
    
    Energies[i] = E_now
    Pressures[i] = P
end

# --- PLOTTING RESULTS ---
p2 = plot(R_scan, Energies, lw=2, color=:blue, legend=false,
          xlabel="Radius R (a.u.)", ylabel="Energy (Ha)", title="Energy vs Confinement")
hline!(p2, [-0.5], color=:black, ls=:dash, label="Free Atom") # Free Hydrogen Limit

p3 = plot(R_scan, Pressures, lw=2, color=:red, legend=false,
          xlabel="Radius R (a.u.)", ylabel="Pressure (a.u.)", title="Pressure on Wall")

# Combine into one figure
l = @layout [a; [b c]]
final_plot = plot(p1, p2, p3, layout=l, size=(900, 850))

display(final_plot)
savefig("confined_hydrogen_pressure.pdf")
