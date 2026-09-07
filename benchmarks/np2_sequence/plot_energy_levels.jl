using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))

using AtomicSplines
using WignerSymbols
using HalfIntegers
using LinearAlgebra
using JLD2
using Plots

const term_3P = LSTerm(1, 2)
const term_1D = LSTerm(2, 0)
const term_1S = LSTerm(0, 0)
const parent_2P = LSTerm(1, 1)

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
        return 0.0
    end
end

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

function compute_p2_SO_reduced_matrix_element(bra_term::LSTerm, ket_term::LSTerm)
    N = 2 
    l = 1; s = HalfInt(1/2)
    L = bra_term.L; two_S = bra_term.two_S; S = two_S // 2
    Lp = ket_term.L; two_Sp = ket_term.two_S; Sp = two_Sp // 2

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
                L_bar = parent_b.L; two_S_bar = parent_b.two_S; S_bar = two_S_bar // 2
                six_j_L = wigner6j(L, 1, Lp, l, L_bar, l)
                six_j_S = wigner6j(S, 1, Sp, s, S_bar, s)
                phase_power = Float64(L_bar) + Float64(S_bar) + l + Float64(s) + L + Float64(Sp)
                phase = iseven(Int(round(phase_power))) ? 1.0 : -1.0
                total_reduced_element += amp_b * amp_k * phase * six_j_L * six_j_S
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
            bra = terms[i]; ket = terms[j]
            reduced_matrix_elem = compute_p2_SO_reduced_matrix_element(bra, ket)
            bra_S = bra.two_S // 2; ket_S = ket.two_S // 2
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

function get_levels(filepath::String, n_val::Int, F2_val::Float64)
    archive = jldopen(filepath, "r")
    E_total = archive["E_total"]
    R_grid = archive["R_grid"]
    V_eff = archive["V_eff"]
    P_np = archive["P_$(n_val)p"]
    close(archive)
    
    zeta_np = compute_zeta(R_grid, V_eff, P_np)
    # The user re-ran ROHF as Configuration Average, so no shift needed
    E_avg = E_total
    
    e0, _ = assemble_and_diagonalize_J_block(0, [term_3P, term_1S], E_avg, F2_val, zeta_np)
    e1, _ = assemble_and_diagonalize_J_block(1, [term_3P], E_avg, F2_val, zeta_np)
    e2, _ = assemble_and_diagonalize_J_block(2, [term_3P, term_1D], E_avg, F2_val, zeta_np)
    
    return [e0..., e1..., e2...]
end

c_cm = [0.0, 29429.5, 20.8, 62.2, 11958.3]
si_cm = [0.0, 17048.7, 113.6, 327.1, 8335.5]
ge_cm = [0.0, 20265.2, 514.5, 1345.5, 9483.4]
sn_cm = [0.0, 17989.0, 1617.1, 3529.6, 11160.6]

labels = ["³P_0", "¹S_0", "³P_1", "³P_2", "¹D_2"]

p = plot(layout=(1,4), size=(1300, 500), title=["Carbono (Z=6)" "Silicio (Z=14)" "Germanio (Z=32)" "Estaño (Z=50)"], legend=false, ylabel="Energía (cm⁻¹)")

for (i, (lvls, name)) in enumerate(zip([c_cm, si_cm, ge_cm, sn_cm], ["Carbono", "Silicio", "Germanio", "Estaño"]))
    sorted_idx = sortperm(lvls)
    sorted_lvls = lvls[sorted_idx]
    sorted_labels = labels[sorted_idx]
    
    x_positions = [0.15, 0.35, 0.55, 0.75, 0.9]
    for (j, lvl) in enumerate(sorted_lvls)
        plot!(p[i], [0, 1], [lvl, lvl], lw=2, color=:blue)
        annotate!(p[i], x_positions[j], lvl, text(sorted_labels[j], 10, :bottom, :black))
    end
    # Set y limits based on the highest level
    ylims!(p[i], -1000, maximum(sorted_lvls) * 1.1)
    xlims!(p[i], -0.2, 1.2)
    xticks!(p[i], :none)
end

savefig("np2_energy_levels.pdf")
println("Saved plot to np2_energy_levels.pdf")
