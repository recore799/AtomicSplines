using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using AtomicSplines, JLD2, LinearAlgebra
include("../../src/ci.jl")

function interaction_coefficient(l1::Int, l2::Int, L::Int, k::Int)
    w6j = wigner6j(Float64, l1, l1, L, l2, l2, k)
    w3j = wigner3j(Float64, l1, k, l2, 0, 0, 0)
    if abs(w3j) < 1e-10 || abs(w6j) < 1e-10
        return 0.0
    end
    rme2 = (2*l1 + 1) * (2*l2 + 1) * w3j^2
    return (-1)^L * w6j * rme2
end

function get_h_core(ws::SolverWorkspace, a::Orbital, b::Orbital, core_orbs::Vector{Orbital})
    if a.l != b.l
        return 0.0
    end
    H_core_mat = ws.T .+ ws.V .+ ((a.l * (a.l + 1)) / 2.0) .* ws.R_inv2
    val = dot(a.coeffs, H_core_mat * b.coeffs)
    for c in core_orbs
        val += 2.0 * (2*c.l + 1) * get_cached_Rk!(ws, a, c, b, c, 0)
        for k in abs(a.l - c.l):(a.l + c.l)
            w3j = wigner3j(Float64, a.l, k, c.l, 0, 0, 0)
            if abs(w3j) > 1e-10
                coeff = (2*c.l + 1) * w3j^2
                val -= coeff * get_cached_Rk!(ws, a, c, c, b, k)
            end
        end
    end
    return val
end

data = load("germanium_rohf_results_R30.0.jld2")
R_max = data["R_max"]
orbitals = data["orbitals"]
num_splines = data["num_splines"]
Z = 32.0
ws = cached_init_scf_workspace(R_max, num_splines, Val(8), Z; γ=3.0)
val_orb = orbitals[end]
core_orbs = orbitals[1:end-1]
build_total_J_matrix!(ws, core_orbs)
J_core = copy(ws.J)
H_core_s = ws.T .+ ws.V .+ ((0.0) / 2.0) .* ws.R_inv2 .+ J_core
assemble_K_matrix!(ws, ws.K_mats[0], 0, core_orbs)
F_s = H_core_s .- ws.K_mats[0]
evals_s, evecs_s = eigen(Symmetric(F_s), Symmetric(ws.S))

active_s = extract_virtuals(evals_s, evecs_s, 0, count(o->o.l==0, core_orbs), 2, 1:ws.basis.num_splines, ws.basis.num_splines, ws)
println("5s h_core = ", get_h_core(ws, active_s[1], active_s[1], core_orbs))
