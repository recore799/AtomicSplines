using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines

using JLD2

function calculate_carbon_term_splitting(filename::String)
    # 1. Load the converged SCF data from the JLD2 file
    data = jldopen(filename, "r")
    orbitals = data["orbitals"]
    E_av = data["E_total"]
    
    # 2. Extract the 2p orbital (it is the 3rd orbital in your array)
    orb_2p = orbitals[3]
    
    # You need to rebuild the workspace to use compute_Rk
    # (Assuming you have a function to quickly reconstruct ws from R_max and N_elems)
    R_max = data["R_max"]
    Z = 6.0
    N_elems = 100
    basis = generate_basis(R_max, N_elems, Val(7), γ=2.5)
    ws = init_scf_workspace(basis, Z)
    
    # 3. Compute F^2(pp) using your R^k function
    # F^2(2p, 2p) = R^2(2p, 2p, 2p, 2p)
    F2_pp = compute_Rk(ws, orb_2p, orb_2p, orb_2p, orb_2p, 2)
    
    println("Configuration Average Energy (E_av) : ", E_av, " Ha")
    println("Quadrupole Slater Integral F^2(pp)  : ", F2_pp, " Ha")
    
    # 4. Calculate the specific LS Term energies (Cowan Eq 12.59)
    E_1S = E_av + (6.0 / 25.0) * F2_pp
    E_1D = E_av + (1.0 / 25.0) * F2_pp
    E_3P = E_av - (3.0 / 25.0) * F2_pp
    
    println("\n--- LS Term Energies ---")
    println("E(^1S) : ", E_1S, " Ha")
    println("E(^1D) : ", E_1D, " Ha")
    println("E(^3P) : ", E_3P, " Ha")
    
    close(data)
    
    return F2_pp, E_1S, E_1D, E_3P
end


calculate_carbon_term_splitting("carbon_rohf_results_R30.jld2")
