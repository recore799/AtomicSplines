using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using JLD2
using Plots

function plot_np2_sequence()
    # 1. Load Data
    carbon_data = jldopen("carbon_rohf_results_R30.0.jld2", "r")
    silicon_data = jldopen("silicon_rohf_results_R30.0.jld2", "r")
    germanium_data = jldopen("germanium_rohf_results_R30.0.jld2", "r")
    
    # We use Carbon's grid as the reference since they are all identical
    r_grid_dense = carbon_data["R_grid"]
    
    # 2. Recreate the Fischer grid limits and transform to sqrt(r)
    # The textbooks typically plot against sqrt(r) to expand the region near the nucleus
    r_ticks = [0.0, 0.125, 0.5, 1.125, 2.0, 3.125, 4.5, 6.125, 8.0, 10.125, 12.5]
    sqrt_r_ticks = sqrt.(r_ticks)
    tick_labels = ["0.0", "0.12", "0.5", "1.12", "2.0", "3.12", "4.5", "6.12", "8.0", "10.1", "12.5"]
    
    sqrt_r_grid = sqrt.(r_grid_dense)
    
    # Filter the grid so we don't plot out to R=30.0 where it's all zero
    max_plot_r = 12.5
    valid_idx = findall(r -> r <= max_plot_r, r_grid_dense)
    
    sqrt_r_plot = sqrt_r_grid[valid_idx]
    
    # --- PLOT 1: VALENCE WAVEFUNCTIONS P_np(r) ---
    p1 = plot(title="Radial Wavefunctions of Valence np² Sequence", 
             xlabel="r (a.u.)", 
             ylabel="P_np(r)", 
             legend=:topright, 
             size=(850, 500),
             margin=5*Plots.mm,
             xticks=(sqrt_r_ticks, tick_labels),
             xlims=(0, sqrt(max_plot_r)))
             
    # Extract wavefunctions and enforce positive phase convention
    P_2p = carbon_data["P_2p"][valid_idx]
    if P_2p[findfirst(x -> abs(x) > 1e-4, P_2p)] < 0; P_2p .*= -1.0; end
    
    P_3p = silicon_data["P_3p"][valid_idx]
    if P_3p[findfirst(x -> abs(x) > 1e-4, P_3p)] < 0; P_3p .*= -1.0; end
    
    P_4p = germanium_data["P_4p"][valid_idx]
    if P_4p[findfirst(x -> abs(x) > 1e-4, P_4p)] < 0; P_4p .*= -1.0; end
    
    plot!(p1, sqrt_r_plot, P_2p, label="Carbon 2p", lw=2.5, color=:blue, linestyle=:solid)
    plot!(p1, sqrt_r_plot, P_3p, label="Silicon 3p", lw=2.5, color=:green, linestyle=:dash)
    plot!(p1, sqrt_r_plot, P_4p, label="Germanium 4p", lw=2.5, color=:red, linestyle=:dot)
    
    display(p1)
    savefig(p1, "np2_valence_wavefunctions.pdf")
    
    # --- PLOT 2: EFFECTIVE POTENTIALS V_eff(r) ---
    p2 = plot(title="Effective Central Potentials for np² Sequence", 
             xlabel="r (a.u.)", 
             ylabel="V_eff(r) (Ha)", 
             legend=:bottomright, 
             size=(850, 500),
             margin=5*Plots.mm,
             xticks=(sqrt_r_ticks, tick_labels),
             xlims=(0, sqrt(max_plot_r)),
             ylims=(-50.0, 5.0)) # Limit y-axis to see the valence region clearly
             
    V_eff_C = carbon_data["V_eff"][valid_idx]
    V_eff_Si = silicon_data["V_eff"][valid_idx]
    V_eff_Ge = germanium_data["V_eff"][valid_idx]
    
    plot!(p2, sqrt_r_plot, V_eff_C, label="Carbon V_eff (Z=6)", lw=2, color=:blue)
    plot!(p2, sqrt_r_plot, V_eff_Si, label="Silicon V_eff (Z=14)", lw=2, color=:green)
    plot!(p2, sqrt_r_plot, V_eff_Ge, label="Germanium V_eff (Z=32)", lw=2, color=:red)
    
    display(p2)
    savefig(p2, "np2_effective_potentials.pdf")
    
    # Close archives
    close(carbon_data)
    close(silicon_data)
    close(germanium_data)
    
    println("Saved np2_valence_wavefunctions.pdf and np2_effective_potentials.pdf")
end

plot_np2_sequence()
