using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using LinearAlgebra
using Printf

function solve_helium_scf(R_box)
    helium = Atom(2.0, [Orbital(1, 0, 1.7)]) 
    orb1s = helium.orbitals[1]

    n_elems = max(80, Int(round(15 * R_box))) 
    k_order = 5
    basis = generate_basis(R_box, n_elems, k_order, γ=2.0)

    T, V_nuc, S = assemble_core(basis, helium.Z)

    solve_orbital!(orb1s, helium, basis)

    MAX_ITER = 60
    MIXING = 0.6
    E_old = 0.0
    E_total = 0.0
    
    for iter in 1:MAX_ITER
        y_coeffs = solve_poisson_potential(basis, orb1s.coeffs, T)
        J_mat = assemble_J_matrix(basis, y_coeffs)

        old_coeffs = copy(orb1s.coeffs)

        solve_orbital!(orb1s, helium, basis, J_mat)

        E_J = dot(orb1s.coeffs, J_mat * orb1s.coeffs)
        E_total = 2 * orb1s.energy - E_J

        orb1s.coeffs = MIXING * orb1s.coeffs + (1.0 - MIXING) * old_coeffs
        n_norm = sqrt(dot(orb1s.coeffs, S * orb1s.coeffs))
        orb1s.coeffs ./= n_norm
        
        diff = abs(E_total - E_old)
        
        if iter > 1 && diff < 1e-9
             break 
        end
        
        E_old = E_total
    end

    return E_total
end


println("\n=======================================================")
println(" VALORES CALCULADOS (HELIO CONFINADO - HF/SCF)")
println("=======================================================\n")

radii_he = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]

println("--- DATOS PARA TABLA HELIO ---")
println("r0\t\tE_0 (Hartree-Fock)")
println("-" ^ 40)

for r in radii_he
    E = solve_helium_scf(r)
    
    println(@sprintf("%.1f\t\t%.4f", r, E))
end
println("\n")
