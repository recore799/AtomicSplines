using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using LinearAlgebra
using Printf

"""
    solve_helium()

Solves the ground state of the Helium atom (Z=2) using the Hartree-Fock method 
expanded in a B-spline basis. Demonstrates the core functionality of AtomicSplines 
including basis generation, matrix assembly, and the Self-Consistent Field (SCF) loop.
"""
function solve_helium()
    println("================================================================")
    println("   HELIUM SOLVER (Hartree-Fock / B-Splines)")
    println("================================================================")

    # 1. System Configuration
    Z = 2.0             # Nuclear charge
    R_MAX = 10.0        # Box size (a.u.)
    N_ELEMS = 60        # Number of finite elements
    ORDER = 7           # B-spline order (k)
    GAMMA = 3.0         # Exponential knot distribution factor
    
    basis = generate_basis(R_MAX, N_ELEMS, ORDER, γ=GAMMA)
    println("Basis: $(basis.num_splines) splines, Order $(basis.order)")

    # 2. Assemble Hamiltonian Operators
    #    T: Kinetic energy, V_nuc: Nuclear potential, S: Overlap
    print("Assembling Core Matrices... ")
    time_core = @elapsed T, V_nuc, S = assemble_core(basis, Z)
    println("Done ($(round(time_core, digits=4)) s)")

    # 3. Initial Guess (He+ Approximation)
    #    Solve H_core * c = E * S * c ignoring electron repulsion
    H_core = T + V_nuc
    
    # Apply Dirichlet Boundary Conditions (u(0)=0, u(R)=0)
    # We solve on the subspace of internal splines [2 : N-1]
    active = 2:(basis.num_splines - 1)
    
    evals, evecs = eigen(H_core[active, active], S[active, active])
    
    # Reconstruct full coefficient vector
    c_current = zeros(Float64, basis.num_splines)
    perm = sortperm(Real.(evals))
    c_current[active] = evecs[:, perm[1]]
    
    # Normalize
    c_current ./= sqrt(dot(c_current, S * c_current))
    
    E_initial = evals[perm[1]]
    println("Initial Guess (He+): $(E_initial) Ha (Target: -2.0)")

    # 4. Self-Consistent Field (SCF) Loop
    MAX_ITER = 30
    MIXING = 0.3      # Linear mixing parameter to stabilize convergence
    TOL = 1e-9        # Energy convergence tolerance
    E_old = 0.0
    
    println("\nStarting SCF Iterations...")
    println("Iter | Total Energy (Ha) | Delta E    | Time (s)")
    println("--------------------------------------------------")
    
    for iter in 1:MAX_ITER
        t_start = time()

        # A. Solve Poisson Equation for Hartree Potential
        #    Returns coefficients for U(r) such that V_H(r) = U(r)/r
        y_coeffs = solve_poisson_potential(basis, c_current, T)
        
        # B. Build Direct Interaction Matrix J
        J = assemble_J_matrix(basis, y_coeffs)
        
        # C. Construct Fock Matrix
        #    F = H_core + J (Restricted Hartree-Fock for closed shell)
        F = H_core + J
        
        # D. Solve Generalized Eigenproblem (F*c = e*S*c)
        evals, evecs = eigen(F[active, active], S[active, active])

        # E. Update Coefficients
        perm = sortperm(Real.(evals))
        c_new = zeros(Float64, basis.num_splines)
        c_new[active] = evecs[:, perm[1]]
        
        # Normalize
        c_new ./= sqrt(dot(c_new, S * c_new))
        
        # F. Compute Total Energy
        #    E_total = 2*epsilon - <J> (double counting correction)
        epsilon = evals[perm[1]]
        E_J = dot(c_new, J * c_new)
        E_total = 2 * epsilon - E_J
        
        # G. Mixing & Convergence Check
        c_current = MIXING * c_current + (1.0 - MIXING) * c_new
        c_current ./= sqrt(dot(c_current, S * c_current)) 
        
        t_iter = time() - t_start
        delta = abs(E_total - E_old)
        
        @printf("%4d | %.10f    | %.2e   | %.4f\n", iter, E_total, delta, t_iter)
        
        if delta < TOL
            println("--------------------------------------------------")
            println("Converged!")
            println("Final Energy: $(E_total) Ha")
            println("Ref Value   : -2.861679995 Ha")
            break
        end
        E_old = E_total
    end
end

solve_helium()
