using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using JLD2
using Printf
using LinearAlgebra

function validate_converged_orbitals(filename::String)
    # 1. Load the strictly optimized ROHF data
    data = jldopen(filename, "r")
    orbitals = data["orbitals"]
    E_total = data["E_total"]
    R_max = data["R_max"]
    Z = 6.0
    
    # 2. Rebuild the B-spline workspace to access the exact analytical matrices
    # We must rebuild because the T (Kinetic) matrix was not saved, and we need it for the Virial Theorem.
    N_elems = 100
    basis = generate_basis(R_max, N_elems, Val(7), γ=2.5)
    ws = init_scf_workspace(basis, Z)
    
    # 3. Compute Total Kinetic Energy (T)
    # The kinetic energy is a one-electron operator. The total kinetic energy of the atom 
    # is the strict sum of the expectation value of the T matrix, weighted by orbital occupancy.
    T_total = 0.0
    for orb in orbitals
        if orb.occ > 0.0
            t_val = dot(orb.coeffs, ws.T * orb.coeffs)
            if orb.l == 1
                t_val += dot(orb.coeffs, ws.V2 * orb.coeffs)
            end
            T_total += orb.occ * t_val
        end
    end
    
    # 4. Compute Total Potential Energy (V) and the Virial Ratio
    # Since E_total = T_total + V_total, we can extract the true effective V.
    V_total = E_total - T_total
    virial_ratio = -V_total / T_total
    
    println("============================================================")
    println("          INTERNAL VALIDATION: THE VIRIAL THEOREM           ")
    println("============================================================")
    @printf(" Total Energy (E)         : %15.8f Ha\n", E_total)
    @printf(" Total Kinetic Energy (T) : %15.8f Ha\n", T_total)
    @printf(" Total Potential (V)      : %15.8f Ha\n", V_total)
    @printf(" Virial Ratio (-V/T)      : %15.8f\n", virial_ratio)
    
    # 5. Compute Radial Expectation Values
    println("\n============================================================")
    println("     RADIAL EXPECTATION VALUES & ORBITAL GEOMETRY           ")
    println("============================================================")
    @printf("%-6s | %-14s | %-14s | %-14s\n", "Orbital", "<r> (Bohr)", "<1/r> (Bohr⁻¹)", "Kinetic <T>")
    println("-"^58)
    
    for orb in orbitals
        if orb.occ > 0.0
            # <r> is extracted directly from your dipole interaction matrix ws.R
            exp_r = dot(orb.coeffs, ws.R * orb.coeffs)
            
            # <1/r> is extracted from the nuclear potential matrix ws.V
            # Since V_ij = <B_i | -Z/r | B_j>, we divide by -Z to isolate <1/r>
            exp_inv_r = dot(orb.coeffs, ws.V * orb.coeffs) / (-Z)
            
            # Individual orbital kinetic energy
            exp_t = dot(orb.coeffs, ws.T * orb.coeffs)
            
            orb_label = string(orb.n, orb.l == 0 ? "s" : (orb.l == 1 ? "p" : "d"))
            @printf("%-6s | %14.8f | %14.8f | %14.8f\n", orb_label, exp_r, exp_inv_r, exp_t)
        end
    end
    println("============================================================")
    
    close(data)
end

# Execute the validation script
validate_converged_orbitals("carbon_rohf_results_R30.jld2")
