using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using JLD2

function calculate_carbon_term_splitting(filename::String)
    # 1. Load the converged SCF data from the JLD2 file
    data = jldopen(filename, "r")
    orbitals = data["orbitals"]
    E_raw_scf = data["E_total"]
    
    # 2. Extract the 2p orbital (it is the 3rd orbital in your array)
    orb_2p = orbitals[3]
    
    # 3. Rebuild the workspace to use compute_Rk
    R_max = data["R_max"]
    Z = 6.0
    N_elems = 100
    basis = generate_basis(R_max, N_elems, Val(7), γ=2.5)
    ws = init_scf_workspace(basis, Z)
    
    # 4. Compute the Monopole (k=0) and Quadrupole (k=2) Slater Integrals
    # F^k(2p, 2p) = R^k(2p, 2p, 2p, 2p)
    F0_pp = compute_Rk(ws, orb_2p, orb_2p, orb_2p, orb_2p, 0)
    F2_pp = compute_Rk(ws, orb_2p, orb_2p, orb_2p, orb_2p, 2)
    
    # 5. Correct the Fractional Occupancy Self-Interaction Error
    # The SCF engine calculated w^2 / 2 = 2.0 pairs.
    # Physics dictates w(w-1) / 2 = 1.0 pair.
    # We must subtract exactly 1.0 spurious pair's worth of energy.
    E_pair = - (3.0 / 25.0) * F2_pp
    
    # The true Configuration Average
    E_av = E_raw_scf - E_pair
    
    println("--- Energy Correction Protocol ---")
    println("Raw SCF Energy (with ghost repulsion) : ", E_raw_scf, " Ha")
    println("Monopole Slater Integral F^0(pp)      : ", F0_pp, " Ha")
    println("Quadrupole Slater Integral F^2(pp)    : ", F2_pp, " Ha")
    println("Spurious Pair Energy Removed          : ", E_pair, " Ha")
    println("Corrected Configuration Average (E_av): ", E_av, " Ha")
    
    # 6. Calculate the specific LS Term energies (Cowan Eq 12.59)
    E_1S = E_av + (6.0 / 25.0) * F2_pp
    E_1D = E_av + (1.0 / 25.0) * F2_pp
    E_3P = E_av - (5.0 / 25.0) * F2_pp
    
    println("\n--- Final LS Term Energies ---")
    println("E(^1S) : ", E_1S, " Ha")
    println("E(^1D) : ", E_1D, " Ha")
    println("E(^3P) : ", E_3P, " Ha")
    
    close(data)
    
    return F0_pp, F2_pp, E_1S, E_1D, E_3P
end

calculate_carbon_term_splitting("carbon_rohf_results_R30.jld2")
