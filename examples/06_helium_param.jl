using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using LinearAlgebra
using Printf

function solve_helium_parametric()
    basis = generate_basis(20.0, 60, Val(7), γ=3.0)
    
    ws = init_scf_workspace(basis, 2.0)

    println("Workspace Initialized. Basis Order K=$(typeof(basis).parameters[1])")

    println("Quadrature Nodes: ", length(basis.gl_nodes))
    println("First Tensor Sum: ", sum(ws.interaction_tensors[div(60, 2)]))

    # Initial Guess
    active = 2:(basis.num_splines - 1)
    H_core = ws.T + ws.V
    evals, evecs = eigen(H_core[active, active], ws.S[active, active])

    # Reconstruct full coefficient vector
    c_current = zeros(Float64, basis.num_splines)
    perm = sortperm(Real.(evals))
    c_current[active] = evecs[:, perm[1]]

    # Normalize
    c_current ./= sqrt(dot(c_current, ws.S * c_current))

    E_initial = evals[perm[1]]
    println("Initial Guess (He+): $(E_initial) Ha (Target: -2.0)")

    # SCF Loop
    MIXING = 0.3
    TOL = 1e-9
    E_old = 0.0

    println("\nStarting SCF Iterations...")
    println("Iter | Total Energy (Ha) | Delta E    | Time (s)")
    println("--------------------------------------------------")
    
    for iter in 1:30
        t_start = time()
        # Solve Poisson (Instant!)
        #    Use the pre-computed Cholesky factorization in ws.poisson_fact
        #    (You'll need to adapt solve_poisson to take the 'ws' struct)
        y_coeffs = solve_poisson_tensor(ws, c_current) 

        
        # Assemble J (Uses ws.scratch_vals, no allocs)
        J = assemble_J_matrix_tensor(ws, y_coeffs)


        # Construct Fock Matrix
        F = H_core + J

        evals, evecs = eigen(F[active, active], ws.S[active, active])

        # Update coefficients
        perm = sortperm(Real.(evals))
        c_new = zeros(Float64, basis.num_splines)
        c_new[active] = evecs[:, perm[1]]

        # Normalize
        c_new ./= sqrt(dot(c_new, ws.S * c_new))

        # Compute total energy
        epsilon = evals[perm[1]]
        E_j = dot(c_new, J * c_new)
        E_total = 2 * epsilon - E_j

        # Mixing and convergence check
        c_current = MIXING * c_current + (1.0 - MIXING) * c_new
        c_current ./= sqrt(dot(c_current, ws.S * c_current))

        t_iter = time() - t_start
        delta = abs(E_total - E_old)

        @printf("%4d | %.10f   | %.2e   | %.4f\n", iter, E_total, delta, t_iter)

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
solve_helium_parametric()
