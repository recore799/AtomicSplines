using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines

using WignerSymbols
using HalfIntegers
using LinearAlgebra
using Printf
using JLD2

# Parent Configuration States (p^2)
const parent_3P = LSTerm(1, 2)
const parent_1D = LSTerm(2, 0)
const parent_1S = LSTerm(0, 0)

# Child Configuration States (p^3)
const child_4S  = LSTerm(0, 3)
const child_2D  = LSTerm(2, 1)
const child_2P  = LSTerm(1, 1)



function compute_zeta(dense_grid::Vector{Float64}, V_eff::Vector{Float64}, P_nl::Vector{Float64})
    N = length(dense_grid)
    dV_dr = zeros(Float64, N)
    integrand = zeros(Float64, N)
    
    # The fine-structure constant (atomic units)
    alpha = 1.0 / 137.035999
    
    # 1. Finite-Difference Numerical Derivative (Second-Order Central)
    dV_dr[1] = (V_eff[2] - V_eff[1]) / (dense_grid[2] - dense_grid[1])
    dV_dr[N] = (V_eff[N] - V_eff[N-1]) / (dense_grid[N] - dense_grid[N-1])
    
    for i in 2:(N-1)
        dV_dr[i] = (V_eff[i+1] - V_eff[i-1]) / (dense_grid[i+1] - dense_grid[i-1])
    end
    
    # 2. Assemble the Spin-Orbit Integrand
    for i in 1:N
        r = dense_grid[i]
        if r > 1e-12
            integrand[i] = (P_nl[i]^2) * (1.0 / r) * dV_dr[i]
        else
            # Analytically enforced limit: P_nl(r)^2 dominates the 1/r^3 divergence for p-orbitals
            integrand[i] = 0.0 
        end
    end
    
    # 3. Trapezoidal Quadrature
    zeta_val = 0.0
    for i in 1:(N-1)
        dr = dense_grid[i+1] - dense_grid[i]
        zeta_val += 0.5 * (integrand[i] + integrand[i+1]) * dr
    end
    
    # Multiply by the relativistic prefactor
    return (alpha^2 / 2.0) * zeta_val
end



"""
    compute_p3_SO_reduced_matrix_element(bra_term::LSTerm, ket_term::LSTerm)

Computes the exact many-body reduced matrix element ⟨p³ LS || V^(11) || p³ L'S'⟩ 
for the Spin-Orbit interaction utilizing Fractional Parentage Coefficients.
"""
function compute_p3_SO_reduced_matrix_element(bra_term::LSTerm, ket_term::LSTerm)
    # The number of equivalent electrons in the active shell
    N = 3 
    
    # Fundamental quantum numbers for a p-electron
    l = 1
    s = HalfInt(1/2)
    
    # Extract the total L and S for the child bra and ket states
    L  = bra_term.L
    two_S = bra_term.two_S
    S = two_S // 2
    Lp = ket_term.L
    two_Sp = ket_term.two_S
    Sp = two_Sp // 2

    # The single-particle reduced matrix elements: ⟨l||l||l⟩ and ⟨s||s||s⟩
    # ⟨l||l||l⟩ = sqrt(l(l+1)(2l+1)) -> For p-electrons (l=1): sqrt(6)
    # ⟨s||s||s⟩ = sqrt(s(s+1)(2s+1)) -> For electrons (s=1/2): sqrt(1.5)
    red_l = sqrt(l * (l + 1.0) * (2l + 1.0))
    red_s = sqrt(Float64(s) * (Float64(s) + 1.0) * (2.0 * Float64(s) + 1.0))
    single_e_product = red_l * red_s

    # The square root dimension prefactor: sqrt((2L+1)(2S+1)(2L'+1)(2S'+1))
    dim_factor = sqrt((2L + 1.0) * (2 * Float64(S) + 1.0) * (2Lp + 1.0) * (2 * Float64(Sp) + 1.0))

    # Retrieve the parent pathways from your established theoretical mapping
    bra_parents = get_parent_amplitudes(bra_term)
    ket_parents = get_parent_amplitudes(ket_term)

    total_reduced_element = 0.0

    # Iterate over the parent states of the bra
    for (parent_b, amp_b) in bra_parents
        # Iterate over the parent states of the ket
        for (parent_k, amp_k) in ket_parents
            
            # The core of the fractional parentage expansion: Orthogonality of the parent states.
            # The matrix element is strictly zero unless the decoupled parent core is identical.
            if parent_b == parent_k
                L_bar = parent_b.L
                two_S_bar = parent_b.two_S
                S_bar = two_S_bar // 2
                
                # 1. Evaluate the Wigner 6-j symbols for Orbital and Spin recoupling
                # { L  1  L' }
                # { l L_bar l }
                six_j_L = wigner6j(L, 1, Lp, l, L_bar, l)
                
                # { S  1  S' }
                # { s S_bar s }
                six_j_S = wigner6j(S, 1, Sp, s, S_bar, s)
                
                # 2. Compute the rigorous parity phase: (-1)^(L_bar + S_bar + l + s + L + S' + 1)
                # We extract the power. Summing half-integers requires care.
                phase_power = Float64(L_bar) + Float64(S_bar) + l + Float64(s) + L + Float64(Sp) + 1.0
                
                # Since phase_power is always an integer for valid physical states,
                # we cast it to Int to safely check parity.
                phase = iseven(Int(round(phase_power))) ? 1.0 : -1.0
                
                # 3. Accumulate the contribution from this specific shared parent pathway
                pathway_contribution = amp_b * amp_k * phase * six_j_L * six_j_S
                
                total_reduced_element += pathway_contribution
            end
        end
    end

    # Multiply by the number of equivalent electrons (N), the dimension factor, and single-particle elements
    final_matrix_element = N * total_reduced_element * dim_factor * single_e_product
    
    return final_matrix_element
end

"""
    assemble_and_diagonalize_J_block(J_target::HalfInt, terms::Vector{LSTerm}, E_avg::Float64, F2::Float64, zeta::Float64)

Constructs the Hamiltonian for a specific J-manifold, applies the electrostatic and 
Spin-Orbit perturbations, and extracts the fine-structure eigenvalues and eigenvectors.
"""
function assemble_and_diagonalize_J_block(J_target::HalfInt, terms::Vector{LSTerm}, E_avg::Float64, F2::Float64, zeta::Float64)
    N = length(terms)
    H_matrix = zeros(Float64, N, N)

    for i in 1:N
        for j in i:N
            bra = terms[i]
            ket = terms[j]
            
            # 1. Spin-Orbit Contribution
            # Calculate the reduced many-body matrix element via Fractional Parentage
            reduced_matrix_elem = compute_p3_SO_reduced_matrix_element(bra, ket)
            
            # Calculate the Wigner 6-j symbol for J-coupling
            two_bra_S = bra.two_S
            bra_S = two_bra_S // 2
            two_ket_S = ket.two_S
            ket_S = two_ket_S // 2

            six_j = wigner6j(bra.L, bra_S, J_target, ket_S, ket.L, 1)
            # Calculate Parity Phase
            phase_val = Float64(bra.L) + Float64(ket_S) + Float64(J_target)
            phase = iseven(Int(round(phase_val))) ? 1.0 : -1.0
            
            # Total SO matrix element
            H_SO = zeta * phase * six_j * reduced_matrix_elem

            if i == j
                # 2. Diagonal Electrostatic Contribution
                c2_coeff = compute_child_F2_coefficient(bra, 3)
                H_elec = E_avg + (c2_coeff * F2)
                
                H_matrix[i, i] = H_elec + H_SO
            else
                # 3. Off-Diagonal Symmetry
                H_matrix[i, j] = H_SO
                H_matrix[j, i] = H_matrix[i, j]
            end
        end
    end

    # Diagonalize the symmetric Hamiltonian
    eigen_decomp = eigen(Symmetric(H_matrix))
    
    return eigen_decomp.values, eigen_decomp.vectors
end

"""
    format_fine_structure_report(E_avg::Float64, F2::Float64, zeta::Float64)

Executes the solver for all J-manifolds and formats the output for publication.
"""
function format_fine_structure_report(E_avg::Float64, F2::Float64, zeta::Float64)
    # Define the allowed states for the J = 3/2 manifold
    terms_J32 = [child_4S, child_2D, child_2P]
    
    println("==================================================")
    println(" Fine Structure Resolution: Phosphorus 3p³")
    println("==================================================")
    
    # Process J = 3/2
    energies_32, vectors_32 = assemble_and_diagonalize_J_block(HalfInt(3/2), terms_J32, E_avg, F2, zeta)
    
    println("\n[ J = 3/2 Manifold ]")
    println("--------------------------------------------------")
    @printf("%-15s | %-20s\n", "Energy (a.u.)", "State Composition (%)")
    println("--------------------------------------------------")
    
    for k in 1:3
        energy = energies_32[k]
        # Calculate percentage composition from the squared eigenvector components
        pct_4S = (vectors_32[1, k]^2) * 100
        pct_2D = (vectors_32[2, k]^2) * 100
        pct_2P = (vectors_32[3, k]^2) * 100
        
        @printf("%.6f        | %5.1f%% ⁴S, %5.1f%% ²D, %5.1f%% ²P\n", energy, pct_4S, pct_2D, pct_2P)
    end
    
    # Process J = 5/2 (Only ²D allowed)
    energies_52, _ = assemble_and_diagonalize_J_block(HalfInt(5/2), [child_2D], E_avg, F2, zeta)
    println("\n[ J = 5/2 Manifold ]")
    @printf("%.6f        | 100.0%% ²D\n", energies_52[1])

    # Process J = 1/2 (Only ²P allowed)
    energies_12, _ = assemble_and_diagonalize_J_block(HalfInt(1/2), [child_2P], E_avg, F2, zeta)
    println("\n[ J = 1/2 Manifold ]")
    @printf("%.6f        | 100.0%% ²P\n", energies_12[1])
    println("==================================================")
end


"""
    execute_phosphorus_fine_structure(filepath::String)

Main execution block: reads the converged SCF data, computes the theoretical 
Slater shifts, integrates the SO parameter, and formats the final report.
"""
function execute_phosphorus_fine_structure(filepath::String)
    println("Loading SCF Configuration Average data from: $filepath")
    
    # 1. Load the variables from the modified JLD2 archive
    archive = jldopen(filepath, "r")
    E_total   = archive["E_total"]
    R_grid    = archive["R_grid"]
    V_eff     = archive["V_eff"]
    P_3p      = archive["P_3p"]
    close(archive)
    
    # 2. Compute Radial Parameter
    zeta_3p = compute_zeta(R_grid, V_eff, P_3p)
    @printf("Calculated Radial Spin-Orbit Parameter (ζ_3p): %.8f Ha\n", zeta_3p)
    
    # 3. Extract the F^2 Integral
    # For a rigorous report, this should be explicitly calculated from the B-splines
    # (e.g., F2 = compute_slater_integral(ws, P_3p, P_3p, 2)). 
    # Assuming a placeholder evaluated value for the execution flow:
    F2_3p = 0.1550  # Replace with actual AtomicSplines F^2 output
    
    # 4. Trigger the previously written diagonalization formatter
    # This invokes assemble_and_diagonalize_J_block internally
    format_fine_structure_report(E_total, F2_3p, zeta_3p)
end


execute_phosphorus_fine_structure("phosphorus_rohf_CA_results_R30.0.jld2")
