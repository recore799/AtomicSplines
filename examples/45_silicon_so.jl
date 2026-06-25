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

# Parent Configuration State for p^1
const parent_2P = LSTerm(1, 1)

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

function get_p2_parent_amplitudes(target_term::LSTerm)
    return [(parent_2P, 1.0)]
end

function get_p2_F2_coefficient(term::LSTerm)
    if term == term_1S
        return 10.0 / 25.0
    elseif term == term_1D
        return 1.0 / 25.0
    elseif term == term_3P
        return -5.0 / 25.0
    else
        error("Unknown p^2 term")
    end
end

function compute_p2_SO_reduced_matrix_element(bra_term::LSTerm, ket_term::LSTerm)
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

    bra_parents = get_p2_parent_amplitudes(bra_term)
    ket_parents = get_p2_parent_amplitudes(ket_term)

    total_reduced_element = 0.0

    for (parent_b, amp_b) in bra_parents
        for (parent_k, amp_k) in ket_parents
            if parent_b == parent_k
                L_bar = parent_b.L
                two_S_bar = parent_b.two_S
                S_bar = two_S_bar // 2
                
                six_j_L = wigner6j(L, 1, Lp, l, L_bar, l)
                six_j_S = wigner6j(S, 1, Sp, s, S_bar, s)
                
                phase_power = Float64(L_bar) + Float64(S_bar) + l + Float64(s) + L + Float64(Sp)
                phase = iseven(Int(round(phase_power))) ? 1.0 : -1.0
                
                pathway_contribution = amp_b * amp_k * phase * six_j_L * six_j_S
                total_reduced_element += pathway_contribution
            end
        end
    end

    return N * total_reduced_element * dim_factor * single_e_product
end

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
                c2_coeff = get_p2_F2_coefficient(bra)
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

function execute_silicon_spin_orbit()
    filepath = "silicon_rohf_results_R30.0.jld2"
    println("Loading SCF Configuration Average data from: $filepath")
    
    archive = jldopen(filepath, "r")
    E_total   = archive["E_total"]
    R_grid    = archive["R_grid"]
    V_eff     = archive["V_eff"]
    P_3p      = archive["P_3p"]
    Z = 14.0
    close(archive)
    
    # Extract Exact F2 parameter calculated in the HF-t state from 41_silicon.jl
    F2_val = 0.165948
    
    # Compute the theoretical Zeta(3p) directly from the V_eff potential
    zeta_computed_ha = compute_zeta(R_grid, V_eff, P_3p)
    ha_to_cm = 219474.6313702
    zeta_computed_cm = zeta_computed_ha * ha_to_cm
    
    # Optional override for exact Froese Fischer value if needed:
    # zeta_np = 149.0 / ha_to_cm # uncomment and change to exact FF value if desired
    zeta_np = zeta_computed_ha
    zeta_cm = zeta_np * ha_to_cm
    
    println("\n==================================================")
    println(" Fine Structure Resolution: Silicon 3p² System")
    println("==================================================")
    @printf("  Term-Dependent Energy (³P)   : %.6f Ha\n", E_total)
    @printf("  Slater F²(3p,3p) Integral    : %.6f Ha\n", F2_val)
    @printf("  Calculated ζ(3p)             : %.3f cm^-1 (%.8f Ha)\n", zeta_cm, zeta_np)
    println("==================================================")

    terms_J0 = [term_3P, term_1S]
    terms_J1 = [term_3P]
    terms_J2 = [term_3P, term_1D]
    
    # The term-dependent SCF gives us the ^3P energy, not the Configuration Average energy!
    # Because E_total in silicon_rohf_results_R30.0.jld2 is the E(³P) energy, we need to shift it
    # to be the E_avg, or we just adjust the F2 coefficients relative to ³P instead of relative to E_avg.
    # In the matrix elements:
    # E(³P) = E_avg - (5/25)*F2  => E_avg = E(³P) + (5/25)*F2
    E_avg = E_total + (5.0 / 25.0) * F2_val
    
    energies_0, vectors_0 = assemble_and_diagonalize_J_block(0, terms_J0, E_avg, F2_val, zeta_np)
    println("\n[ J = 0 Manifold ]")
    println("--------------------------------------------------")
    @printf("%-15s | %-20s\n", "Energy (a.u.)", "State Composition (%)")
    println("--------------------------------------------------")
    for k in 1:2
        pct_3P = (vectors_0[1, k]^2) * 100
        pct_1S = (vectors_0[2, k]^2) * 100
        @printf("%.6f        | %5.1f%% ³P, %5.1f%% ¹S\n", energies_0[k], pct_3P, pct_1S)
    end
    
    energies_1, _ = assemble_and_diagonalize_J_block(1, terms_J1, E_avg, F2_val, zeta_np)
    println("\n[ J = 1 Manifold ]")
    println("--------------------------------------------------")
    @printf("%.6f        | 100.0%% ³P\n", energies_1[1])

    energies_2, vectors_2 = assemble_and_diagonalize_J_block(2, terms_J2, E_avg, F2_val, zeta_np)
    println("\n[ J = 2 Manifold ]")
    println("--------------------------------------------------")
    for k in 1:2
        pct_3P = (vectors_2[1, k]^2) * 100
        pct_1D = (vectors_2[2, k]^2) * 100
        @printf("%.6f        | %5.1f%% ³P, %5.1f%% ¹D\n", energies_2[k], pct_3P, pct_1D)
    end
    
    println("==================================================")
    
    E_3P_0 = energies_0[1] * ha_to_cm
    E_3P_1 = energies_1[1] * ha_to_cm
    E_3P_2 = energies_2[1] * ha_to_cm
    
    split_1_0 = E_3P_1 - E_3P_0
    split_2_1 = E_3P_2 - E_3P_1
    
    println("\n--- Theoretical Splitting of the ³P Multiplet ---")
    @printf("  ³P_1 - ³P_0 : %7.3f cm^-1 (Expected ~%.3f cm^-1)\n", split_1_0, zeta_cm / 2.0)
    @printf("  ³P_2 - ³P_1 : %7.3f cm^-1 (Expected ~%.3f cm^-1)\n", split_2_1, zeta_cm)
    println("==================================================")
end

execute_silicon_spin_orbit()
