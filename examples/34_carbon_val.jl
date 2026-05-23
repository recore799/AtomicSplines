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

    ws = cached_init_scf_workspace(R_max, N_elems, Val(7), Z; γ=2.5, calc_R_matrices=true)
    basis = ws.basis

    # 3. Compute Total Kinetic Energy (T)
    # The kinetic energy is a one-electron operator. The total kinetic energy of the atom 
    # is the strict sum of the expectation value of the T matrix, weighted by orbital occupancy.
    T_total = 0.0
    for orb in orbitals
        if orb.occ > 0.0
            t_val = dot(orb.coeffs, ws.T * orb.coeffs)
            if orb.l == 1
                t_val += dot(orb.coeffs, ws.R_inv2 * orb.coeffs)
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
    @printf("%-4s | %-12s | %-12s | %-12s | %-12s\n", 
            "Orb", "<r> (Bohr)", "<r²> (Bohr²)", "<1/r> (Bohr⁻¹)", "<1/r³> (Bohr⁻³)")
    println("-"^74)
    
    for orb in orbitals
        if orb.occ > 0.0
            # <r> is extracted directly from your dipole interaction matrix ws.R
            exp_r = dot(orb.coeffs, ws.R * orb.coeffs)

            exp_r2 = dot(orb.coeffs, ws.R2 * orb.coeffs)
            
            # <1/r> is extracted from the nuclear potential matrix ws.V
            # Since V_ij = <B_i | -Z/r | B_j>, we divide by -Z to isolate <1/r>
            exp_inv_r = dot(orb.coeffs, ws.V * orb.coeffs) / (-Z)

            exp_inv_r3 = dot(orb.coeffs, ws.R_inv3 * orb.coeffs)

            # Individual orbital kinetic energy
            exp_t = dot(orb.coeffs, ws.T * orb.coeffs)
            
            orb_label = string(orb.n, orb.l == 0 ? "s" : (orb.l == 1 ? "p" : "d"))
            @printf("%-4s | %12.6f | %12.6f | %12.6f | %12.6f\n", 
                    orb_label, exp_r, exp_r2, exp_inv_r, exp_inv_r3)
        end
    end
    println("==========================================================================")

    # 6. Extract the Valence 4p Orbital and Compute Slater Integrals
    # The 4p orbital strictly resides at index 8 in the Germanium configuration
    orb_2p = orbitals[3]
    
    F0_2p = compute_Rk(ws, orb_2p, orb_2p, orb_2p, orb_2p, 0)
    F2_2p = compute_Rk(ws, orb_2p, orb_2p, orb_2p, orb_2p, 2)
    
    println("\n==========================================================================")
    println("              VALENCE 2p SLATER INTEGRALS (TERM-DEPENDENT)                ")
    println("==========================================================================")
    @printf(" Monopole Interaction   F⁰(2p, 2p) : %10.8f Ha\n", F0_2p)
    @printf(" Quadrupole Interaction F²(2p, 2p) : %10.8f Ha\n", F2_2p)
    println("==========================================================================")
 
    
    close(data)
end

# Execute the validation script
validate_converged_orbitals("carbon_rohf_results_R30.jld2")
