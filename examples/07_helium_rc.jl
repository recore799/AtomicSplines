using Pkg
Pkg.activate(joinpath(@__DIR__, "..")) 

using AtomicSplines
using LinearAlgebra
using Printf
using Roots

function solve_helium_scf(R_box)
    helium = Atom(2.0, [Orbital(1, 0, 1.7)]) 
    orb1s = helium.orbitals[1]

    n_elems = max(40, Int(round(15 * R_box))) 
    k_order = 5
    basis = generate_basis(R_box, n_elems, k_order, γ=2.0)

    T, V_nuc, S = assemble_core(basis, helium.Z)

    # Initial Guess
    solve_orbital!(orb1s, helium, basis)

    MAX_ITER = 60
    MIXING = 0.6
    E_old = 0.0
    E_total = 0.0
    
    for iter in 1:MAX_ITER
        y_coeffs = solve_poisson_potential(basis, orb1s.coeffs, T)
        J_mat = assemble_J_matrix(basis, y_coeffs)

        old_coeffs = copy(orb1s.coeffs)

        # Solve eigenvalue problem (h + J)c = εSc
        solve_orbital!(orb1s, helium, basis, J_mat)

        # Calculate Total Energy (Standard HF)
        E_J = dot(orb1s.coeffs, J_mat * orb1s.coeffs)
        E_total = 2 * orb1s.energy - E_J

        # Mix coefficients (Linear mixing)
        orb1s.coeffs = MIXING * orb1s.coeffs + (1.0 - MIXING) * old_coeffs
        
        # Re-normalize
        n_norm = sqrt(dot(orb1s.coeffs, S * orb1s.coeffs))
        orb1s.coeffs ./= n_norm
        
        diff = abs(E_total - E_old)
        if iter > 1 && diff < 1e-8
             break 
        end
        E_old = E_total
    end

    return E_total
end

function find_critical_radius(r_guess_min, r_guess_max)
    # Define function f(r) such that f(r) = 0 at critical radius
    f(r) = solve_helium_scf(r)
    
    println("Searching for critical radius between $r_guess_min and $r_guess_max...")
    
    return find_zero(f, (r_guess_min, r_guess_max), Roots.Bisection())
end

println("\n=======================================================")
println(" CRITICAL RADIUS CALCULATION (Helium - 1s Orbital)")
println("=======================================================\n")

# For Orbital Energy = 0, rc is usually between 0.5 and 1.0 for He.
# If seeking Total Energy = 0, try range (0.1, 0.6)
rc = find_critical_radius(0.9, 1.2)

println("Critical Radius (rc):\t" * @sprintf("%.8f a.u.", rc))
