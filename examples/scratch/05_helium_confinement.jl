using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using LinearAlgebra
using Printf

"""
    solve_helium_confined(R_box)

Runs the Helium SCF calculation for a specific hard-wall box radius.
Returns the converged Total Energy (Ha).
"""
function solve_helium_confined(R_box)
    # 1. System Configuration
    Z = 2.0            
    N_ELEMS = 60       
    ORDER = 7          
    # Slightly lower Gamma for confinement to ensure we have resolution 
    # near the wall (where the wavefunction is forced to zero).
    GAMMA = 2.0        
    
    basis = generate_basis(R_box, N_ELEMS, Val(ORDER), γ=GAMMA)

    # 2. Assemble Hamiltonian Operators
    #    Silence the print output for the loop
    ws = init_scf_workspace(basis, Z)

    # 3. Initial Guess (He+)
    H_core = ws.T + ws.V
    active = 2:(basis.num_splines - 1) # Hard Wall BCs
    
    evals, evecs = eigen(H_core[active, active], ws.S[active, active])
    
    c_current = zeros(Float64, basis.num_splines)
    perm = sortperm(Real.(evals))
    c_current[active] = evecs[:, perm[1]]
    c_current ./= sqrt(dot(c_current, ws.S * c_current))
    
    # 4. SCF Loop
    MAX_ITER = 50 # Increased slightly for hard compression
    MIXING = 0.4  # Conservative mixing
    TOL = 1e-9    
    E_old = 0.0
    E_final = 0.0
    
    for iter in 1:MAX_ITER
        # A. Poisson
        y_coeffs = solve_poisson_J(ws, c_current)
        
        # B. J Matrix
        J = assemble_J_matrix(ws, y_coeffs)
        
        # C. Fock
        F = H_core + J
        
        # D. Solve
        evals, evecs = eigen(F[active, active], ws.S[active, active])

        # E. Update
        perm = sortperm(Real.(evals))
        c_new = zeros(Float64, basis.num_splines)
        c_new[active] = evecs[:, perm[1]]
        c_new ./= sqrt(dot(c_new, ws.S * c_new))
        
        # F. Energy
        epsilon = evals[perm[1]]
        E_J = dot(c_new, J * c_new)
        E_total = 2 * epsilon - E_J
        
        # G. Mixing
        c_current = MIXING * c_current + (1.0 - MIXING) * c_new
        c_current ./= sqrt(dot(c_current, ws.S * c_current)) 
        
        delta = abs(E_total - E_old)
        E_final = E_total
        
        if delta < TOL
            return E_final
        end
        E_old = E_total
    end
    
    println("Warning: R=$R_box did not fully converge (Delta=$delta)")
    return E_final
end

# ==============================================================================
# MAIN EXECUTION LOOP
# ==============================================================================

radii = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]

println("\n==================================================================")
println(" HELIUM CONFINEMENT SCAN (Hard Wall)")
println("==================================================================")
println(" Radius (a.u.) | Total Energy (Ha) | Delta vs Free Limit")
println("------------------------------------------------------------------")

# Free Helium Limit (Infinite Box)
E_free = -2.861679995 

for R in radii
    # Run the solver
    E = solve_helium_confined(R)
    
    diff = E - E_free
    
    # Format the output nicely
    @printf(" %4.1f          | %12.8f      | +%.6f\n", R, E, diff)
end
println("------------------------------------------------------------------")
