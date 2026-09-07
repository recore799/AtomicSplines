using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using JLD2
using Plots
using AtomicSplines

function get_evaluated_orbs(data_file, K, gamma)
    data = jldopen(data_file, "r")
    r_grid = data["R_grid"]
    max_plot_r = 12.5
    valid_idx = findall(r -> r <= max_plot_r, r_grid)
    sqrt_r_plot = sqrt.(r_grid[valid_idx])
    
    N_splines = data["num_splines"]
    N_elems = N_splines - K
    R_max = data["R_max"]
    basis = AtomicSplines.generate_basis(R_max, N_elems, Val(K); γ=gamma)
    
    orbs = data["orbitals"]
    evals = Dict{String, Vector{Float64}}()
    for orb in orbs
        l_str = ['s', 'p', 'd', 'f'][orb.l + 1]
        name = "$(orb.n)$(l_str)"
        P_r = evaluate_orbital(basis, orb.coeffs, r_grid)
        P_r = P_r[valid_idx]
        if P_r[findfirst(x -> abs(x) > 1e-4, P_r)] < 0
            P_r .*= -1.0
        end
        evals[name] = P_r
    end
    close(data)
    return evals, sqrt_r_plot
end

function plot_all_orbitals()
    c_evals, sqrt_r_plot = get_evaluated_orbs(joinpath(@__DIR__, "carbon_rohf_results_R30.0.jld2"), 7, 2.50)
    si_evals, _ = get_evaluated_orbs(joinpath(@__DIR__, "silicon_rohf_results_R30.0.jld2"), 7, 2.50)
    ge_evals, _ = get_evaluated_orbs(joinpath(@__DIR__, "germanium_rohf_results_R30.0.jld2"), 8, 3.00)
    sn_evals, _ = get_evaluated_orbs(joinpath(@__DIR__, "..", "..", "benchmarks", "np2_sequence", "tin_rohf_results_3P_R30.0.jld2"), 8, 4.00)
    
    max_plot_r = 12.5
    r_ticks = [0.0, 0.125, 0.5, 1.125, 2.0, 3.125, 4.5, 6.125, 8.0, 10.125, 12.5]
    sqrt_r_ticks = sqrt.(r_ticks)
    tick_labels = ["0.0", "0.12", "0.5", "1.12", "2.0", "3.12", "4.5", "6.12", "8.0", "10.1", "12.5"]
    
    p_c = plot(title="Carbono (Z=6)", ylabel="P_nl(r)", margin=3*Plots.mm,
               xticks=(sqrt_r_ticks, tick_labels), xlims=(0, sqrt(max_plot_r)), legend=:topright)
    for (name, P_r) in sort(collect(c_evals))
        plot!(p_c, sqrt_r_plot, P_r, label=name, lw=2)
    end
    
    p_si = plot(title="Silicio (Z=14)", ylabel="P_nl(r)", margin=3*Plots.mm,
               xticks=(sqrt_r_ticks, tick_labels), xlims=(0, sqrt(max_plot_r)), legend=:topright)
    for (name, P_r) in sort(collect(si_evals))
        plot!(p_si, sqrt_r_plot, P_r, label=name, lw=2)
    end
    
    p_ge = plot(title="Germanio (Z=32)", ylabel="P_nl(r)", margin=3*Plots.mm,
               xticks=(sqrt_r_ticks, tick_labels), xlims=(0, sqrt(max_plot_r)), legend=:topright)
    for (name, P_r) in sort(collect(ge_evals))
        plot!(p_ge, sqrt_r_plot, P_r, label=name, lw=2)
    end
    
    p_sn = plot(title="Estaño (Z=50)", xlabel="r (u.a.)", ylabel="P_nl(r)", margin=3*Plots.mm,
               xticks=(sqrt_r_ticks, tick_labels), xlims=(0, sqrt(max_plot_r)), legend=:topright)
    for (name, P_r) in sort(collect(sn_evals))
        plot!(p_sn, sqrt_r_plot, P_r, label=name, lw=2)
    end
    
    p_all = plot(p_c, p_si, p_ge, p_sn, layout=(4,1), size=(850, 1300))
    savefig(p_all, joinpath(@__DIR__, "..", "..", "docs", "figures", "np2_all_orbitals.pdf"))
    println("Saved docs/figures/np2_all_orbitals.pdf")
end

plot_all_orbitals()
