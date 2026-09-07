using Pkg
Pkg.activate(joinpath(@__DIR__, "..")) 


using AtomicSplines, LinearAlgebra, Printf

function tune_grid_neon()
    println("--- NEON BASIS RESOLUTION TEST (Target: -50.00 Ha) ---")
    println("N_splines | R_max | Gamma |  E_1s (Hydrogenic) |  Diff")
    println("-"^60)

    # Test a few configurations
    configs = [
        (50, 10.0, 3.0), # Old settings (maybe sparse?)
        (120, 8.0,  3.0), # Your recent failed attempt
        (120, 8.0,  4.0), # Sharper focus near nucleus
        (150, 10.0, 4.0), # High density
        (200, 15.0, 5.0)  # Overkill (Should be perfect)
    ]

    for (N, R, G) in configs
        basis = generate_basis(R, N, 7, γ=G)
        
        # Assemble ONLY the one-electron Hamiltonian for Z=10
        # We need higher quadrature here too!
        # (Assuming standard assemble_core uses order+2, which might be low for Z=10)
        T, V_nuc, S = assemble_core(basis, 10.0)
        
        # Solve Eigenvalues
        vals, _ = eigen(Symmetric(Matrix((T+V_nuc)[2:end-1, 2:end-1])), 
                        Symmetric(Matrix(S[2:end-1, 2:end-1])))
        
        e_calc = vals[1]
        diff = abs(e_calc - (-50.0))
        
        @printf("%9d | %5.1f | %5.1f | %18.10f | %.1e", N, R, G, e_calc, diff)
        if diff < 0.1
            println("  [PASS]")
        else
            println("  [FAIL]")
        end
    end
end

tune_grid_neon()
