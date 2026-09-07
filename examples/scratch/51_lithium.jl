using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using JLD2
using Printf
using LinearAlgebra

function validate_lithium_orbitals(filename::String)
    # 1. Load the strictly optimized HF-t data
    data = jldopen(filename, "r")
    orbitals = data["orbitals"]
    E_total = data["E_total"]
    R_max = data["R_max"]
    Z = 14.0
    
    # 2. Rebuild the B-spline workspace
    # Assumes init_scf_workspace now also provides ws.R2 and ws.R_inv3
    N_elems = 100
    basis = generate_basis(R_max, N_elems, Val(7), γ=2.5)
    ws = init_scf_workspace(basis, Z)
    
    # 3. Compute Total Kinetic Energy (T)
    # CRITICAL FIX: The centrifugal term (ws.V2) must be added for all l > 0 states.
    T_total = 0.0
    for orb in orbitals
        if orb.occ > 0.0
            if orb.l == 0
                t_val = dot(orb.coeffs, ws.T * orb.coeffs)
            else
                t_val = dot(orb.coeffs, (ws.T + ws.V2) * orb.coeffs)
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
    @printf(" Total Energy (E)         : %15.8f Ha\n", E_total)
    @printf(" Total Kinetic Energy (T) : %15.8f Ha\n", T_total)
    @printf(" Total Potential (V)      : %15.8f Ha\n", V_total)
    @printf(" Virial Ratio (-V/T)      : %15.8f\n", virial_ratio)
    
    # 5. Compute Radial Expectation Values
    println("\n==========================================================================")
    println("               RADIAL EXPECTATION VALUES & ORBITAL GEOMETRY               ")
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
    
    # 6. Extract the Valence 3p Orbital and Compute Slater Integrals
    # orb_3p = orbitals[5]
    
    # F0_3p = compute_Rk(ws, orb_3p, orb_3p, orb_3p, orb_3p, 0)
    # F2_3p = compute_Rk(ws, orb_3p, orb_3p, orb_3p, orb_3p, 2)
    
    # println("\n==========================================================================")
    # println("              VALENCE 3p SLATER INTEGRALS (TERM-DEPENDENT)                ")
    # println("==========================================================================")
    # @printf(" Monopole Interaction   F⁰(3p, 3p) : %10.6f Ha\n", F0_3p)
    # @printf(" Quadrupole Interaction F²(3p, 3p) : %10.6f Ha\n", F2_3p)
    # println("==========================================================================")
    
    close(data)
end

# Execute the validation script
validate_lithium_orbitals("lithium_rohf_results_R30.0.jld2")
