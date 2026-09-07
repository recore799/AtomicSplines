using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using JLD2
using Plots

function plot_np2_sequence()
    # 1. Load Data
    carbon_data = jldopen(joinpath(@__DIR__, "carbon_rohf_results_R30.0.jld2"), "r")
    silicon_data = jldopen(joinpath(@__DIR__, "silicon_rohf_results_R30.0.jld2"), "r")
    germanium_data = jldopen(joinpath(@__DIR__, "germanium_rohf_results_R30.0.jld2"), "r")
    tin_data = jldopen(joinpath(@__DIR__, "..", "..", "benchmarks", "np2_sequence", "tin_rohf_results_3P_R30.0.jld2"), "r")
    
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
    p1 = plot(title="Funciones de Onda Radiales de Valencia (Secuencia np²)", 
             xlabel="r (u.a.)", 
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
    
    P_5p = tin_data["P_5p"][valid_idx]
    if P_5p[findfirst(x -> abs(x) > 1e-4, P_5p)] < 0; P_5p .*= -1.0; end
    
    plot!(p1, sqrt_r_plot, P_2p, label="Carbono 2p", lw=2.5, color=:blue, linestyle=:solid)
    plot!(p1, sqrt_r_plot, P_3p, label="Silicio 3p", lw=2.5, color=:green, linestyle=:dash)
    plot!(p1, sqrt_r_plot, P_4p, label="Germanio 4p", lw=2.5, color=:red, linestyle=:dot)
    plot!(p1, sqrt_r_plot, P_5p, label="Estaño 5p", lw=2.5, color=:orange, linestyle=:dashdot)
    
    display(p1)
    savefig(p1, joinpath(@__DIR__, "..", "..", "docs", "figures", "np2_valence_wavefunctions.pdf"))
    
    # --- PLOT 2: EFFECTIVE POTENTIALS V_eff(r) (Pseudo-log) ---
    pseudo_log(V) = sign(V) * log10(1 + abs(V)/0.1)
    y_ticks_vals = [0, -0.5, -1, -5, -10, -50, -100, -500, -1000, -3500]
    y_ticks_pos = pseudo_log.(y_ticks_vals)
    y_tick_labels = string.(y_ticks_vals)
    
    p2 = plot(title="Potenciales Centrales Efectivos (Secuencia np²)", 
             xlabel="r (u.a.)", 
             ylabel="V_eff(r) (Ha)", 
             legend=:bottomright, 
             size=(850, 500),
             margin=5*Plots.mm,
             xticks=(sqrt_r_ticks, tick_labels),
             yticks=(y_ticks_pos, y_tick_labels),
             xlims=(0, sqrt(max_plot_r)),
             ylims=(pseudo_log(-3600), pseudo_log(2)))
             
    V_eff_C = carbon_data["V_eff"][valid_idx]
    V_eff_Si = silicon_data["V_eff"][valid_idx]
    V_eff_Ge = germanium_data["V_eff"][valid_idx]
    V_eff_Sn = tin_data["V_eff"][valid_idx]
    
    plot!(p2, sqrt_r_plot, pseudo_log.(V_eff_C), label="Carbono V_eff (Z=6)", lw=2, color=:blue)
    plot!(p2, sqrt_r_plot, pseudo_log.(V_eff_Si), label="Silicio V_eff (Z=14)", lw=2, color=:green)
    plot!(p2, sqrt_r_plot, pseudo_log.(V_eff_Ge), label="Germanio V_eff (Z=32)", lw=2, color=:red)
    plot!(p2, sqrt_r_plot, pseudo_log.(V_eff_Sn), label="Estaño V_eff (Z=50)", lw=2, color=:orange)
    
    display(p2)
    savefig(p2, joinpath(@__DIR__, "..", "..", "docs", "figures", "np2_effective_potentials.pdf"))

    # --- PLOT 3: RADIAL EFFECTIVE POTENTIALS V_rad(r) = V_eff(r) + 1/r² ---
    r_vals = carbon_data["R_grid"][valid_idx]
    
    # Avoid division by zero at r=0. The first point is usually r=0 in the grid.
    centrifugal = zeros(length(r_vals))
    for i in 1:length(r_vals)
        if r_vals[i] > 1e-10
            centrifugal[i] = 1.0 / (r_vals[i]^2) # l=1, l(l+1)/2 = 1
        else
            centrifugal[i] = 1e6 # cap to avoid infinity in plot
        end
    end
    
    p3 = plot(title="Potenciales Radiales Efectivos (V_eff + 1/r²)", 
             xlabel="r (u.a.)", 
             ylabel="V_rad(r) (Ha)", 
             legend=:bottomright, 
             size=(850, 500),
             margin=5*Plots.mm,
             xticks=(sqrt_r_ticks, tick_labels),
             xlims=(0, sqrt(max_plot_r)),
             ylims=(-300, 50))
             
    V_rad_C = V_eff_C .+ centrifugal
    V_rad_Si = V_eff_Si .+ centrifugal
    V_rad_Ge = V_eff_Ge .+ centrifugal
    V_rad_Sn = V_eff_Sn .+ centrifugal
    
    plot!(p3, sqrt_r_plot, V_rad_C, label="Carbono V_rad (Z=6)", lw=2, color=:blue)
    plot!(p3, sqrt_r_plot, V_rad_Si, label="Silicio V_rad (Z=14)", lw=2, color=:green)
    plot!(p3, sqrt_r_plot, V_rad_Ge, label="Germanio V_rad (Z=32)", lw=2, color=:red)
    plot!(p3, sqrt_r_plot, V_rad_Sn, label="Estaño V_rad (Z=50)", lw=2, color=:orange)
    
    # Horizontal line at E=0
    hline!(p3, [0], color=:black, lw=1, linestyle=:dash, label="")
    
    display(p3)
    savefig(p3, joinpath(@__DIR__, "..", "..", "docs", "figures", "np2_radial_potentials.pdf"))
    
    # Close archives
    close(carbon_data)
    close(silicon_data)
    close(germanium_data)
    close(tin_data)
    
    println("Saved np2_valence_wavefunctions.pdf and np2_effective_potentials.pdf")
end

plot_np2_sequence()
