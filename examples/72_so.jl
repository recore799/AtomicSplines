using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using WignerSymbols
using HalfIntegers
using LinearAlgebra
using Printf
using JLD2

# Allowed Configuration States for p^2
const term_3P = LSTerm(1, 2)
const term_1D = LSTerm(2, 0)
const term_1S = LSTerm(0, 0)

"""
    compute_zeta(dense_grid::Vector{Float64}, V_eff::Vector{Float64}, P_nl::Vector{Float64})

Computes the radial spin-orbit parameter using second-order finite differences 
and trapezoidal quadrature.
"""
function compute_zeta(dense_grid::Vector{Float64}, V_eff::Vector{Float64}, P_nl::Vector{Float64})
    N = length(dense_grid)
    dV_dr = zeros(Float64, N)
    integrand = zeros(Float64, N)
    
    alpha = 1.0 / 137.035999
    
    dV_dr[1] = (V_eff[2] - V_eff[1]) / (dense_grid[2] - dense_grid[1])
    dV_dr[N] = (V_eff[N] - V_eff[N-1]) / (dense_grid[N] - dense_grid[N-1])
    
    for i in 2:(N-1)
        dV_dr[i] = (V_eff[i+1] - V_eff[i-1]) / (dense_grid[i+1] - dense_grid[i-1])
    end
    
    for i in 1:N
        r = dense_grid[i]
        if r > 1e-12
            integrand[i] = (P_nl[i]^2) * (1.0 / r) * dV_dr[i]
        else
            integrand[i] = 0.0 
        end
    end
    
    zeta_val = 0.0
    for i in 1:(N-1)
        dr = dense_grid[i+1] - dense_grid[i]
        zeta_val += 0.5 * (integrand[i] + integrand[i+1]) * dr
    end
    
    return (alpha^2 / 2.0) * zeta_val
end

"""
    compute_p2_SO_reduced_matrix_element(bra_term::LSTerm, ket_term::LSTerm)

Computes the exact many-body reduced matrix element for two equivalent p-electrons.
"""
function compute_p2_SO_reduced_matrix_element(bra_term::LSTerm, ket_term::LSTerm)
    # The number of equivalent electrons is strictly 2 for C, Si, Ge
    N = 2 
    
    l = 1
    s = HalfInt(1/2)
    
    L  = bra_term.L
    two_S = bra_term.two_S
    S = two_S // 2
    Lp = ket_term.L
    two_Sp = ket_term.two_S
    Sp = two_Sp // 2

    red_l = sqrt(l * (l + 1.0) * (2l + 1.0))
    red_s = sqrt(Float64(s) * (Float64(s) + 1.0) * (2.0 * Float64(s) + 1.0))
    single_e_product = red_l * red_s

    dim_factor = sqrt((2L + 1.0) * (2 * Float64(S) + 1.0) * (2Lp + 1.0) * (2 * Float64(Sp) + 1.0))

    bra_parents = get_parent_amplitudes(bra_term)
    ket_parents = get_parent_amplitudes(ket_term)

    total_reduced_element = 0.0

    for (parent_b, amp_b) in bra_parents
        for (parent_k, amp_k) in ket_parents
            if parent_b == parent_k
                L_bar = parent_b.L
                two_S_bar = parent_b.two_S
                S_bar = two_S_bar // 2
                
                six_j_L = wigner6j(L, 1, Lp, l, L_bar, l)
                six_j_S = wigner6j(S, 1, Sp, s, S_bar, s)
                
                phase_power = Float64(L_bar) + Float64(S_bar) + l + Float64(s) + L + Float64(Sp) + 1.0
                phase = iseven(Int(round(phase_power))) ? 1.0 : -1.0
                
                pathway_contribution = amp_b * amp_k * phase * six_j_L * six_j_S
                total_reduced_element += pathway_contribution
            end
        end
    end

    return N * total_reduced_element * dim_factor * single_e_product
end

"""
    assemble_and_diagonalize_J_block(J_target::Int, terms::Vector{LSTerm}, E_avg::Float64, F2::Float64, zeta::Float64)

Constructs the Hamiltonian for a specific integer J-manifold and returns eigenvalues and eigenvectors.
"""
function assemble_and_diagonalize_J_block(J_target::Int, terms::Vector{LSTerm}, E_avg::Float64, F2::Float64, zeta::Float64)
    N = length(terms)
    H_matrix = zeros(Float64, N, N)

    for i in 1:N
        for j in i:N
            bra = terms[i]
            ket = terms[j]
            
            reduced_matrix_elem = compute_p2_SO_reduced_matrix_element(bra, ket)
            
            bra_S = bra.two_S // 2
            ket_S = ket.two_S // 2

            six_j = wigner6j(bra.L, bra_S, J_target, ket_S, ket.L, 1)
            
            phase_val = Float64(bra.L) + Float64(ket_S) + Float64(J_target)
            phase = iseven(Int(round(phase_val))) ? 1.0 : -1.0
            
            H_SO = zeta * phase * six_j * reduced_matrix_elem

            if i == j
                # Utilize 2 electrons for the F2 coefficient calculation
                c2_coeff = compute_child_F2_coefficient(bra, 2)
                H_elec = E_avg + (c2_coeff * F2)
                H_matrix[i, i] = H_elec + H_SO
            else
                H_matrix[i, j] = H_SO
                H_matrix[j, i] = H_matrix[i, j]
            end
        end
    end

    eigen_decomp = eigen(Symmetric(H_matrix))
    return eigen_decomp.values, eigen_decomp.vectors
end

"""
    format_fine_structure_report(element::String, E_avg::Float64, F2::Float64, zeta::Float64)

Executes the solver for the integer J-manifolds associated with a p^2 system.
"""
function format_fine_structure_report(element::String, E_avg::Float64, F2::Float64, zeta::Float64)
    # Define the interacting states within each allowed J-manifold for p^2
    terms_J0 = [term_3P, term_1S]
    terms_J1 = [term_3P]
    terms_J2 = [term_3P, term_1D]
    
    println("==================================================")
    println(" Fine Structure Resolution: $element np² System")
    println("==================================================")
    
    # Process J = 0 Manifold (³P_0 and ¹S_0 mix)
    energies_0, vectors_0 = assemble_and_diagonalize_J_block(0, terms_J0, E_avg, F2, zeta)
    println("\n[ J = 0 Manifold ]")
    println("--------------------------------------------------")
    @printf("%-15s | %-20s\n", "Energy (a.u.)", "State Composition (%)")
    println("--------------------------------------------------")
    for k in 1:2
        pct_3P = (vectors_0[1, k]^2) * 100
        pct_1S = (vectors_0[2, k]^2) * 100
        @printf("%.6f        | %5.1f%% ³P, %5.1f%% ¹S\n", energies_0[k], pct_3P, pct_1S)
    end
    
    # Process J = 1 Manifold (Only ³P_1 allowed)
    energies_1, _ = assemble_and_diagonalize_J_block(1, terms_J1, E_avg, F2, zeta)
    println("\n[ J = 1 Manifold ]")
    println("--------------------------------------------------")
    @printf("%.6f        | 100.0%% ³P\n", energies_1[1])

    # Process J = 2 Manifold (³P_2 and ¹D_2 mix)
    energies_2, vectors_2 = assemble_and_diagonalize_J_block(2, terms_J2, E_avg, F2, zeta)
    println("\n[ J = 2 Manifold ]")
    println("--------------------------------------------------")
    for k in 1:2
        pct_3P = (vectors_2[1, k]^2) * 100
        pct_1D = (vectors_2[2, k]^2) * 100
        @printf("%.6f        | %5.1f%% ³P, %5.1f%% ¹D\n", energies_2[k], pct_3P, pct_1D)
    end
    
    println("==================================================")
    
    # Optional post-processing analytical fraction evaluation for exact state splitting
    exact_splitting_factor = 2.0 / 3.0
    approx_splitting = (energies_2[1] - energies_0[1]) * exact_splitting_factor
    @printf("Rigorous Splitting Evaluation (J=2 - J=0) * 2/3: %.8f Ha\n", approx_splitting)
    println("==================================================")
end

"""
    execute_np2_fine_structure(element::String, filepath::String, F2_val::Float64)

Main execution block for generalized p^2 elements.
"""
function execute_np2_fine_structure(element::String, filepath::String, F2_val::Float64)
    println("Loading SCF Configuration Average data from: $filepath")
    
    # Establish the principal quantum number mapped to the element
    element_n_map = Dict("C" => 2, "Si" => 3, "Ge" => 4)
    if !haskey(element_n_map, element)
        error("Element not supported in this mapping. Please use C, Si, or Ge.")
    end
    n_val = element_n_map[element]
    orbital_key = "P_$(n_val)p"
    
    archive = jldopen(filepath, "r")
    E_total   = archive["E_total"]
    R_grid    = archive["R_grid"]
    V_eff     = archive["V_eff"]
    P_np      = archive[orbital_key]
    close(archive)
    
    zeta_np = compute_zeta(R_grid, V_eff, P_np)
    @printf("Calculated Radial Spin-Orbit Parameter (ζ_%dp): %.8f Ha\n", n_val, zeta_np)
    
    format_fine_structure_report(element, E_total, F2_val, zeta_np)
end

# Example Execution
execute_np2_fine_structure("C", "carbon_rohf_results_R30.0.jld2", 0.2200)
execute_np2_fine_structure("Si", "silicon_rohf_results_R30.0.jld2", 0.1450)
execute_np2_fine_structure("Ge", "germanium_rohf_results_R30.0.jld2", 0.162724)

