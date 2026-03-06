using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using LinearAlgebra
using JLD2
using Plots

# This will load 'orbitals', 'E_total', and 'R_max' back into your workspace
@load "neon_results_R20.0.jld2" orbitals E_total R_max



function enforce_positive_phase!(orbitals::Vector{Orbital})
    for orb in orbitals
        # The first non-zero B-spline coefficient near the origin
        first_active_idx = orb.l + 2 
        
        # If it starts negative, flip the entire wavefunction
        if orb.coeffs[first_active_idx] < 0.0
            orb.coeffs .*= -1.0
        end
    end
end


function plot_orbitals(basis::BSplineBasis, orbitals::Vector{Orbital}, R_max::Float64)
    # Recreate the Fischer grid limits
    # The textbook stops around r = 8.0 to 9.0, so we limit our viewing window.
    r_ticks = [0.0, 0.125, 0.5, 1.125, 2.0, 3.125, 4.5, 6.125, 8.0]
    sqrt_r_ticks = sqrt.(r_ticks)
    
    # Format labels to match the book (e.g., "1.12" instead of "1.125")
    tick_labels = ["0.0", "0.12", "0.5", "1.12", "2.0", "3.12", "4.5", "6.12", "8.0"]
    
    # Define a very dense fine grid uniformly spaced in sqrt(r)
    # This ensures smooth curves when we plot against sqrt(r)
    sqrt_r_points = collect(range(0.0, sqrt(8.5), length=1500))
    r_points = sqrt_r_points .^ 2  # The actual physical r values to evaluate
    
    # Match the exact order of the orbitals array
    # labels = ["1s", "2s", "3s", "4s", "5s", "6s", "2p", "3p", "4p", "5p", "6p", "3d", "4d", "5d" , "4f"]
    labels = ["1s", "2s", "2p"]
    
    # Initialize the plot with the custom X-axis ticks
    p = plot(title="Hartree-Fock Radial Wavefunctions of Neon", 
             xlabel="r (a.u.)", 
             ylabel="P(r)", 
             legend=:outertopright, 
             size=(850, 600),
             margin=5*Plots.mm,
             xticks=(sqrt_r_ticks, tick_labels),
             xlims=(0, sqrt(8.5)))
             
    # Extract and plot each orbital
    for (i, orb) in enumerate(orbitals)
        # Evaluate using the true 'r_points'
        psi_vals = evaluate_orbital(basis, orb.coeffs, r_points)
        
        # Use different line styles based on angular momentum
        ls = if orb.l == 0
            :solid
        elseif orb.l == 1
            :dash
        elseif orb.l == 2
            :dot
        else
            :dashdot
        end
        
        # The key trick: we plot against `sqrt_r_points` instead of `r_points`!
        plot!(p, sqrt_r_points, psi_vals, label=labels[i], lw=2, linestyle=ls)
    end
    
    display(p)
    savefig(p, "neon_orbitals.pdf")
end

function enforce_positive_phase!(orbitals::Vector{Orbital})
    for orb in orbitals
        # The first non-zero B-spline coefficient near the origin
        first_active_idx = orb.l + 2 
        
        # If it starts negative, flip the entire wavefunction
        if orb.coeffs[first_active_idx] < 0.0
            orb.coeffs .*= -1.0
        end
    end
end

function plot_density(basis::BSplineBasis, orbitals::Vector{Orbital}, R_max::Float64)
    # 1. Recreate the Fischer grid limits to match your previous plot
    # r_ticks = [0.0, 0.125, 0.5, 1.125, 2.0, 3.125, 4.5, 6.125, 8.0]
    r_ticks = [0.0, 0.5, 2.0, 4.5, 8.0]
    sqrt_r_ticks = sqrt.(r_ticks)
    tick_labels = ["0.0", "0.12", "0.5", "1.12", "2.0", "3.12", "4.5", "6.12", "8.0"]
    
    # Dense grid mapped to sqrt(r)
    sqrt_r_points = collect(range(0.0, sqrt(8.5), length=1500))
    r_points = sqrt_r_points .^ 2 
    
    # 2. Initialize an array to hold the total density
    total_density = zeros(Float64, length(r_points))
    
    # 3. Accumulate D(r) from each orbital
    for orb in orbitals
        # Evaluate P_nl(r)
        psi_vals = evaluate_orbital(basis, orb.coeffs, r_points)
        # 3. Accumulate D(r) from each orbital

        total_density .+= orb.occ .* (psi_vals .^ 2) 

    end
    
    # 4. Create the plot
    p = plot(title="Total Radial Electron Density of Neon", 
             xlabel="r (a.u.)", 
             ylabel="D(r)", 
             legend=false, 
             size=(850, 500),
             margin=5*Plots.mm,
             xticks=(sqrt_r_ticks, tick_labels),
             xlims=(0, sqrt(8.5)))
             
    # Plot with a solid line and a subtle fill under the curve
    plot!(p, sqrt_r_points, total_density, lw=2, color=:black, fill=(0, 0.2, :blue))
    
    display(p)
    
    savefig(p, "neon_total_density.pdf")
end

N_elems = 100
Z = 36.0
R_max = 20.0

basis = generate_basis(R_max, N_elems, Val(7), γ=2.5)

enforce_positive_phase!(orbitals)
plot_orbitals(basis, orbitals, R_max)
plot_density(basis, orbitals, R_max)
