using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using WignerSymbols
using HalfIntegers
using LinearAlgebra
using Printf

const term_3P = (L=1, two_S=2)
const term_1D = (L=2, two_S=0)
const term_1S = (L=0, two_S=0)

function compute_p2_SO_reduced_matrix_element(bra_term, ket_term)
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

    L_bar = 1
    S_bar = 1//2
    six_j_L = wigner6j(L, 1, Lp, l, L_bar, l)
    six_j_S = wigner6j(S, 1, Sp, s, S_bar, s)
    
    phase_power = Float64(L_bar) + Float64(S_bar) + l + Float64(s) + L + Float64(Sp)
    phase = iseven(Int(round(phase_power))) ? 1.0 : -1.0
    
    pathway_contribution = 1.0 * 1.0 * phase * six_j_L * six_j_S
    total_reduced_element = pathway_contribution

    return N * total_reduced_element * dim_factor * single_e_product
end

function get_H_SO(J_target, bra, ket, zeta)
    reduced_matrix_elem = compute_p2_SO_reduced_matrix_element(bra, ket)
    bra_S = bra.two_S // 2
    ket_S = ket.two_S // 2
    six_j = wigner6j(bra.L, bra_S, J_target, ket_S, ket.L, 1)
    phase_val = Float64(bra.L) + Float64(ket_S) + Float64(J_target)
    phase = iseven(Int(round(phase_val))) ? 1.0 : -1.0
    return zeta * phase * six_j * reduced_matrix_elem
end

println("J=0: <3P0|H|3P0> = ", get_H_SO(0, term_3P, term_3P, 1.0), " zeta")
println("J=0: <3P0|H|1S0> = ", get_H_SO(0, term_3P, term_1S, 1.0), " zeta")
println("J=1: <3P1|H|3P1> = ", get_H_SO(1, term_3P, term_3P, 1.0), " zeta")
println("J=2: <3P2|H|3P2> = ", get_H_SO(2, term_3P, term_3P, 1.0), " zeta")
println("J=2: <3P2|H|1D2> = ", get_H_SO(2, term_3P, term_1D, 1.0), " zeta")
