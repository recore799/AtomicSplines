using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using LinearAlgebra
using Printf

function solve_beryllium()
    println("=== Beryllium (Z=4) RHF Solver ===")
    
    # 1. Setup
    R_max = 15.0
    N_elems = 150
    Z = 4.0
    # Use Order 7 (standard for high precision)
    basis = generate_basis(R_max, N_elems, Val(7), γ=2.5)
    ws = init_scf_workspace(basis, Z)
    
    n = basis.num_splines
    active = 2:(n-1)
    
    # 2. Initial Guess (H-core)
    # Solve H c = E S c
    println("Computing Initial Guess...")
    H_core = ws.T + ws.V
    
    # Dense Solve for safety in guess
    evals, evecs = eigen(Matrix(H_core[active, active]), Matrix(ws.S[active, active]))
    
    # Construct initial orbitals (normalized)
    c_1s = zeros(Float64, n); c_1s[active] = evecs[:, 1]
    c_2s = zeros(Float64, n); c_2s[active] = evecs[:, 2]
    
    # Normalize
    c_1s ./= sqrt(dot(c_1s, ws.S * c_1s))
    c_2s ./= sqrt(dot(c_2s, ws.S * c_2s))
    
    # Density Matrix P (for damping)
    P = 2.0 * (c_1s * c_1s') + 2.0 * (c_2s * c_2s')
    
    E_old = 0.0
    MIXING = 0.3
    
    println("Starting SCF Loop...")
    println("Iter | Energy (Ha)    | Delta    | Time")
    println("----------------------------------------")
    
    for iter in 1:40
        t0 = time()
        
        # Rho_total = 2 * |1s|^2 + 2 * |2s|^2 (APPROXIMATION?)
        rho_vec = 2.0 .* (c_1s .^ 2) .+ 2.0 .* (c_2s .^ 2)
       
        y_1s = solve_poisson_J(ws, c_1s)
        y_2s = solve_poisson_J(ws, c_2s)
        

        y_tot = 2.0 .* y_1s .+ 2.0 .* y_2s
        J_mat = assemble_J_matrix(ws, y_tot)
        
        # Exchange Interaction (K)
        # K_tot = K_1s + K_2s
        K_mat = assemble_K_matrix(ws, [c_1s, c_2s])
        
        # Fock Matrix
        F = H_core + J_mat - K_mat
        
        # Diagonalize
        evals_f, evecs_f = eigen(Matrix(F[active, active]), Matrix(ws.S[active, active]))
        
        # Update Coefficients
        c_1s_new = zeros(Float64, n); c_1s_new[active] = evecs_f[:, 1]
        c_2s_new = zeros(Float64, n); c_2s_new[active] = evecs_f[:, 2]
        
        # Normalize
        c_1s_new ./= sqrt(dot(c_1s_new, ws.S * c_1s_new))
        c_2s_new ./= sqrt(dot(c_2s_new, ws.S * c_2s_new))
        
        h_11 = dot(c_1s_new, H_core * c_1s_new)
        h_22 = dot(c_2s_new, H_core * c_2s_new)
        eps_1 = evals_f[1]
        eps_2 = evals_f[2]
        
        # E_total = (eps_1 + h_11) * 2/2 ? No.
        # RHF Energy: E = sum_i (epsilon_i + h_ii)
        E_total = (eps_1 + h_11) + (eps_2 + h_22)
        
        delta = abs(E_total - E_old)
        t_iter = time() - t0
        
        @printf("%4d | %.8f | %.2e | %.4fs\n", iter, E_total, delta, t_iter)
        
        if delta < 1e-9
            println("Converged!")
            println("Final Energy: $E_total Ha")
            println("Reference   : -14.5730 Ha (approx HF limit)")
            break
        end
        
        # Mixing (Linear Mixing of Coefficients)
        c_1s = MIXING * c_1s + (1-MIXING) * c_1s_new
        c_2s = MIXING * c_2s + (1-MIXING) * c_2s_new
        c_1s ./= sqrt(dot(c_1s, ws.S * c_1s))
        c_2s ./= sqrt(dot(c_2s, ws.S * c_2s))
        
        E_old = E_total
    end
end

solve_beryllium()
