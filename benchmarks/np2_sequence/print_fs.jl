using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using AtomicSplines, JLD2, LinearAlgebra
include("../../src/ci.jl")
data = load("carbon_rohf_results_R30.0.jld2")
orbitals = data["orbitals"]
Z = 6.0
ws = cached_init_scf_workspace(data["R_max"], 100, Val(7), Z; γ=2.5)
core_orbs = orbitals[1:end-1]
build_total_J_matrix!(ws, core_orbs)
J_core = copy(ws.J)
H_core_s = ws.T .+ ws.V .+ J_core
assemble_K_matrix!(ws, ws.K_mats[0], 0, core_orbs)
F_s = H_core_s .- ws.K_mats[0]
evals_s, evecs_s = eigen(Symmetric(F_s), Symmetric(ws.S))
println("F_s evals (first 10):")
println(evals_s[1:10])
