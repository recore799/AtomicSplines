using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using JLD2
using Plots
using AtomicSplines

function get_total_density(data_file, K, gamma)
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
    
    total_density = zeros(Float64, length(valid_idx))
    for orb in orbs
        P_r = evaluate_orbital(basis, orb.coeffs, r_grid)
        P_r = P_r[valid_idx]
        total_density .+= orb.occ .* (P_r .^ 2)
    end
    close(data)
    return total_density, sqrt_r_plot
end

function plot_all_densities()
    c_dens, sqrt_r_plot = get_total_density(joinpath(@__DIR__, "carbon_rohf_results_R30.0.jld2"), 7, 2.50)
    si_dens, _ = get_total_density(joinpath(@__DIR__, "silicon_rohf_results_R30.0.jld2"), 7, 2.50)
    ge_dens, _ = get_total_density(joinpath(@__DIR__, "germanium_rohf_results_R30.0.jld2"), 8, 3.00)
    sn_dens, _ = get_total_density(joinpath(@__DIR__, "..", "..", "benchmarks", "np2_sequence", "tin_rohf_results_3P_R30.0.jld2"), 8, 4.00)
    
    max_plot_r = 12.5
    r_ticks = [0.0, 0.125, 0.5, 1.125, 2.0, 3.125, 4.5, 6.125, 8.0, 10.125, 12.5]
    sqrt_r_ticks = sqrt.(r_ticks)
    tick_labels = ["0.0", "0.12", "0.5", "1.12", "2.0", "3.12", "4.5", "6.12", "8.0", "10.1", "12.5"]
    
    p_c = plot(title="Densidad Total Carbono (Z=6)", ylabel="D(r)", margin=3*Plots.mm,
               xticks=(sqrt_r_ticks, tick_labels), xlims=(0, sqrt(max_plot_r)), legend=false)
    plot!(p_c, sqrt_r_plot, c_dens, lw=2, color=:black, fill=(0, 0.2, :blue))
    
    p_si = plot(title="Densidad Total Silicio (Z=14)", ylabel="D(r)", margin=3*Plots.mm,
               xticks=(sqrt_r_ticks, tick_labels), xlims=(0, sqrt(max_plot_r)), legend=false)
    plot!(p_si, sqrt_r_plot, si_dens, lw=2, color=:black, fill=(0, 0.2, :green))
    
    p_ge = plot(title="Densidad Total Germanio (Z=32)", ylabel="D(r)", margin=3*Plots.mm,
               xticks=(sqrt_r_ticks, tick_labels), xlims=(0, sqrt(max_plot_r)), legend=false)
    plot!(p_ge, sqrt_r_plot, ge_dens, lw=2, color=:black, fill=(0, 0.2, :red))
    
    p_sn = plot(title="Densidad Total Estaño (Z=50)", xlabel="r (u.a.)", ylabel="D(r)", margin=3*Plots.mm,
               xticks=(sqrt_r_ticks, tick_labels), xlims=(0, sqrt(max_plot_r)), legend=false)
    plot!(p_sn, sqrt_r_plot, sn_dens, lw=2, color=:black, fill=(0, 0.2, :orange))
    
    p_all = plot(p_c, p_si, p_ge, p_sn, layout=(4,1), size=(850, 1300))
    savefig(p_all, joinpath(@__DIR__, "..", "..", "docs", "figures", "np2_total_densities.pdf"))
    println("Saved docs/figures/np2_total_densities.pdf")
end

plot_all_densities()
