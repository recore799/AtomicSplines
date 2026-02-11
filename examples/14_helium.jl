using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using LinearAlgebra
using Printf

function solve_helium_old_api_benchmark()
    println("================================================================")
    println("   HELIUM SOLVER (Old API Benchmark)")
    println("================================================================")

    # 1. Setup (Exact same parameters as Optimized version for fair comparison)
    R_MAX = 20.0
    N_ELEMS = 60
    ORDER = 7 
    Z = 2.0
    GAMMA = 3.0
    
    basis = generate_basis(R_MAX, N_ELEMS, ORDER, γ=GAMMA)
    println("Basis: $(basis.num_splines) splines, Order $(basis.order)")

    # 2. Fixed Operators (T, V, S) - Using OLD assemble_core
    print("Assembling Core Matrices... ")
    time_core = @elapsed T, V_nuc, S = assemble_core(basis, Z)
    println("Done ($(round(time_core, digits=4)) s)")

    # 3. Initial Guess (He+ without interaction)
    # Define system using the Old Structs
    helium = Atom(Z, [Orbital(1, 0, Z)]) 
    orb1s = helium.orbitals[1]

    # Solve initial state
    solve_orbital!(orb1s, helium, basis)
    
    # Calculate Initial Energy (2 * epsilon)
    E_initial = orb1s.energy * 2.0
    println("Initial Guess (He+): $(E_initial) Ha")

    # 4. SCF Loop
    MAX_ITER = 30
    MIXING = 0.3 
    TOL = 1e-9
    E_old = 0.0
    
    println("\nStarting SCF (Old API)...")
    println("Iter | Total Energy (Ha) | Delta E    | Time (s)")
    println("--------------------------------------------------")
    
    for iter in 1:MAX_ITER
        t_start = time()

        # A. Poisson (Old API)
        # Note: Old API usually required manual solving of Poisson
        y_coeffs = solve_poisson_potential(basis, orb1s.coeffs, T)
        
        # B. Build J Matrix (Old API)
        J_mat = assemble_J_matrix(basis, y_coeffs)

        # C. Store old coefficients for mixing
        coeffs_prev = copy(orb1s.coeffs)

        # D. Solve Eigenproblem (Old API: solve_orbital!)
        # This function internally builds F = H + J and solves
        solve_orbital!(orb1s, helium, basis, J_mat)

        # E. Calculate Total Energy
        # E_tot = 2*epsilon - <J>
        # Note: orb1s.energy is updated inside solve_orbital!
        E_J = dot(orb1s.coeffs, J_mat * orb1s.coeffs)
        E_total = 2 * orb1s.energy - E_J

        # F. Mixing (Manual implementation to match Optimized logic)
        # Formula: c_current = MIX * c_old + (1-MIX) * c_new
        # In this flow, orb1s.coeffs currently holds c_new
        orb1s.coeffs = (1.0 - MIXING) * orb1s.coeffs + MIXING * coeffs_prev
        
        # Re-normalize (Critical step in manual mixing)
        n_val = sqrt(dot(orb1s.coeffs, S * orb1s.coeffs))
        orb1s.coeffs ./= n_val

        # Timing and Delta
        t_iter = time() - t_start
        delta = abs(E_total - E_old)
        
        @printf("%4d | %.10f    | %.2e   | %.4f\n", iter, E_total, delta, t_iter)
        
        if delta < TOL
            println("--------------------------------------------------")
            println("Converged!")
            println("Final Energy: $(E_total) Ha")
            println("Ref Value   : -2.86168 Ha")
            break
        end
        E_old = E_total
    end
end

solve_helium_old_api_benchmark()
