using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))

using AtomicSplines
using JLD2
using Printf
using LinearAlgebra

function validate_tin_orbitals(filename::String)
    # 1. Load the strictly optimized HF-t data
    data = jldopen(joinpath(@__DIR__, filename), "r")
    orbitals = data["orbitals"]
    E_total = data["E_total"]
    R_max = data["R_max"]
    Z = 50.0
    
    # 2. Rebuild the B-spline workspace
    N_elems = 500
    ws = cached_init_scf_workspace(R_max, N_elems, Val(8), Z; γ=4.0, calc_R_matrices=true)
    
    # 3. Compute Total Kinetic Energy (T)
    T_total = 0.0
    for orb in orbitals
        if orb.occ > 0.0
            if orb.l == 0
                t_val = dot(orb.coeffs, ws.T * orb.coeffs)
            elseif orb.l == 1
                t_val = dot(orb.coeffs, (ws.T + ws.R_inv2) * orb.coeffs)
            elseif orb.l == 2
                t_val = dot(orb.coeffs, (ws.T + 3.0 * ws.R_inv2) * orb.coeffs)
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
            exp_r = dot(orb.coeffs, ws.R * orb.coeffs)
            exp_r2 = dot(orb.coeffs, ws.R2 * orb.coeffs)
            exp_inv_r = dot(orb.coeffs, ws.V * orb.coeffs) / (-Z)
            exp_inv_r3 = dot(orb.coeffs, ws.R_inv3 * orb.coeffs)
            
            orb_label = string(orb.n, orb.l == 0 ? "s" : (orb.l == 1 ? "p" : "d"))
            @printf("%-4s | %12.6f | %12.6f | %12.6f | %12.6f\n", 
                    orb_label, exp_r, exp_r2, exp_inv_r, exp_inv_r3)
        end
    end
    println("==========================================================================")
    
    # 6. Extract the Valence 5p Orbital and Compute Slater Integrals
    # In our tin configuration, 5p is orbital index 11
    orb_5p = orbitals[11]
    
    F0_5p = compute_Rk(ws, orb_5p, orb_5p, orb_5p, orb_5p, 0)
    F2_5p = compute_Rk(ws, orb_5p, orb_5p, orb_5p, orb_5p, 2)
    
    println("\n==========================================================================")
    println("              VALENCE 5p SLATER INTEGRALS (TERM-DEPENDENT)                ")
    println("==========================================================================")
    @printf(" Monopole Interaction   F⁰(5p, 5p) : %10.6f Ha\n", F0_5p)
    @printf(" Quadrupole Interaction F²(5p, 5p) : %10.6f Ha\n", F2_5p)
    println("==========================================================================")
    
    close(data)
end

if !isinteractive()
    if isfile(joinpath(@__DIR__, "tin_rohf_results_3P_R30.0.jld2"))
        validate_tin_orbitals("tin_rohf_results_3P_R30.0.jld2")
    elseif isfile(joinpath(@__DIR__, "tin_rohf_results_av_R30.0.jld2"))
        validate_tin_orbitals("tin_rohf_results_av_R30.0.jld2")
    else
        println("Could not find Tin ROHF results. Please run tin_rohf.jl first.")
    end
end
