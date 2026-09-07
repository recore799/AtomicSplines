using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using LinearAlgebra
using Printf
using JLD2


function validate_germanium_orbitals(filename::String)
    # 1. Load the strictly optimized HF-t data
    data = jldopen(filename, "r")
    orbitals = data["orbitals"]
    E_total = data["E_total"]
    R_max = data["R_max"]
    
    # Germanium specific term and atomic number
    term_symbol = haskey(data, "term") ? data["term"] : "^3P"
    Z = 32.0
    
    # 2. Rebuild the B-spline workspace
    N_elems = 300
    ws = cached_init_scf_workspace(R_max, N_elems, Val(8), Z; γ=3.0, calc_R_matrices=true)
    
    # 3. Compute Total Kinetic Energy (T)
    # RIGOROUS FIX: The centrifugal term l(l+1)/2 must be properly scaled for s, p, and d shells.
    T_total = 0.0
    for orb in orbitals
        if orb.occ > 0.0
            if orb.l == 0
                t_val = dot(orb.coeffs, ws.T * orb.coeffs)
            elseif orb.l == 1
                t_val = dot(orb.coeffs, (ws.T + ws.R_inv2) * orb.coeffs)
            elseif orb.l == 2
                # d-orbital: l=2 -> 2(3)/2 = 3.0
                t_val = dot(orb.coeffs, (ws.T + 3.0 * ws.R_inv2) * orb.coeffs)
            else
                error("Orbital angular momentum l > 2 not supported in this script.")
            end
            T_total += orb.occ * t_val
        end
    end
    
    # 4. Compute Total Potential Energy (V) and the Virial Ratio
    V_total = E_total - T_total
    virial_ratio = -V_total / T_total
    
    println("==========================================================================")
    println("                 INTERNAL VALIDATION: THE VIRIAL THEOREM                  ")
    println("==========================================================================")
    @printf(" Element                  : Germanium (Z = %.1f)\n", Z)
    @printf(" Target State             : 4p^2 (%s)\n", term_symbol)
    @printf(" Total Energy (E)         : %15.8f Ha\n", E_total)
    @printf(" Total Kinetic Energy (T) : %15.8f Ha\n", T_total)
    @printf(" Total Potential (V)      : %15.8f Ha\n", V_total)
    @printf(" Virial Ratio (-V/T)      : %15.8f\n", virial_ratio)
    
    # 5. Compute Radial Expectation Values
    println("\n==========================================================================")
    println("                RADIAL EXPECTATION VALUES & ORBITAL GEOMETRY               ")
    println("==========================================================================")
    @printf("%-4s | %-12s | %-12s | %-12s | %-12s\n", 
            "Orb", "<r> (Bohr)", "<r²> (Bohr²)", "<1/r> (Bohr⁻¹)", "<1/r³> (Bohr⁻³)")
    println("-"^74)
    
    for orb in orbitals
        if orb.occ > 0.0
            # Dipole <r> 
            exp_r = dot(orb.coeffs, ws.R * orb.coeffs)
            
            # Quadrupole <r^2>
            exp_r2 = dot(orb.coeffs, ws.R2 * orb.coeffs)
            
            # Inverse Radius <1/r> (Extracted from the Z-scaled nuclear potential)
            exp_inv_r = dot(orb.coeffs, ws.V * orb.coeffs) / (-Z)
            
            # Inverse Cube <1/r^3> 
            exp_inv_r3 = dot(orb.coeffs, ws.R_inv3 * orb.coeffs)
            
            orb_label = string(orb.n, orb.l == 0 ? "s" : (orb.l == 1 ? "p" : "d"))
            @printf("%-4s | %12.6f | %12.6f | %12.6f | %12.6f\n", 
                    orb_label, exp_r, exp_r2, exp_inv_r, exp_inv_r3)
        end
    end
    println("==========================================================================")
    
    # 6. Extract the Valence 4p Orbital and Compute Slater Integrals
    # The 4p orbital strictly resides at index 8 in the Germanium configuration
    orb_4p = orbitals[8]
    
    F0_4p = compute_Rk(ws, orb_4p, orb_4p, orb_4p, orb_4p, 0)
    F2_4p = compute_Rk(ws, orb_4p, orb_4p, orb_4p, orb_4p, 2)
    
    println("\n==========================================================================")
    println("              VALENCE 4p SLATER INTEGRALS (TERM-DEPENDENT)                ")
    println("==========================================================================")
    @printf(" Monopole Interaction   F⁰(4p, 4p) : %10.8f Ha\n", F0_4p)
    @printf(" Quadrupole Interaction F²(4p, 4p) : %10.8f Ha\n", F2_4p)
    println("==========================================================================")
    
    close(data)
end

validate_germanium_orbitals("germanium_rohf_results_R30.0.jld2")
