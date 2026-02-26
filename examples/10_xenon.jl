using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using LinearAlgebra
using Plots
using Printf
using JLD2


function solve_xenon(R_max; verbose::Bool=true)
    println("=== Xenon (Z=54) ===")
    
    time1 = time()
    N_elems = 120
    Z = 54.0
    basis = generate_basis(R_max, N_elems, Val(7), γ=2.5)
    ws = init_scf_workspace(basis, Z)

    n = basis.num_splines

    # Enforce origin boundary conditions: P(r) ~ r^(l+1)
    active_s = 2:(n-1)  # P(0) = 0
    active_p = 3:(n-1)  # P(0) = 0, P'(0) = 0
    active_d = 4:(n-1)  # P(0) = 0, P'(0) = 0, P''(0) = 0
    
    if verbose
        println("\n--- Información de la Base ---")
        println("  Radio máx (R_max)  : $R_max a.u.")
        println("  Elementos          : $N_elems")
        println("  Splines totales (n): $n")
        println("  Funciones s activas: $(length(active_s))")
        println("  Funciones p activas: $(length(active_p))")
        println("  Funciones d activas: $(length(active_d))")
        println("------------------------------\n")
    end

    # Initialize Orbitals (n, l, occupancy)
    orbitals = [
        Orbital(1, 0, 2.0),  # 1s (idx 1)
        Orbital(2, 0, 2.0),  # 2s (idx 2)
        Orbital(3, 0, 2.0),  # 3s (idx 3)
        Orbital(4, 0, 2.0),  # 4s (idx 4)
        Orbital(5, 0, 2.0),  # 5s (idx 5)
        Orbital(2, 1, 6.0),  # 2p (idx 6)
        Orbital(3, 1, 6.0),  # 3p (idx 7)
        Orbital(4, 1, 6.0),  # 4p (idx 8)
        Orbital(5, 1, 6.0),  # 5p (idx 9)
        Orbital(3, 2, 10.0), # 3d (idx 10)
        Orbital(4, 2, 10.0)  # 4d (idx 11)
    ]
    
    # --- Core Hamiltonians ---
    H_core_s = ws.T + ws.V
    H_core_p = ws.T + ws.V + ws.V2 
    H_core_d = ws.T + ws.V + 3.0 * ws.V2  # l(l+1)/2 = 6/2 = 3
    
    # --- Initial Guess ---
    evals_s, evecs_s = eigen(Symmetric(H_core_s[active_s, active_s]), ws.S[active_s, active_s])
    evals_p, evecs_p = eigen(Symmetric(H_core_p[active_p, active_p]), ws.S[active_p, active_p])
    evals_d, evecs_d = eigen(Symmetric(H_core_d[active_d, active_d]), ws.S[active_d, active_d])
    
    # Assign initial guesses based on principal quantum number
    for (i, orb_idx) in enumerate(1:5) # 1s to 5s
        orbitals[orb_idx].coeffs = zeros(Float64, n)
        orbitals[orb_idx].coeffs[active_s] = evecs_s[:, i]
    end
    for (i, orb_idx) in enumerate(6:9) # 2p to 5p
        orbitals[orb_idx].coeffs = zeros(Float64, n)
        orbitals[orb_idx].coeffs[active_p] = evecs_p[:, i]
    end
    for (i, orb_idx) in enumerate(10:11) # 3d to 4d
        orbitals[orb_idx].coeffs = zeros(Float64, n)
        orbitals[orb_idx].coeffs[active_d] = evecs_d[:, i]
    end
    
    # Normalize initial guesses
    for orb in orbitals
        orb.coeffs ./= sqrt(dot(orb.coeffs, ws.S * orb.coeffs))
    end

    # --- SCF Loop ---
    E_old = 0.0
    MIXING = 0.3
    
    println("Comenzando ciclo SCF...")
    if verbose
        @printf("%-4s | %-14s | %-10s | %-8s\n", 
                "Iter", "E_total (Ha)", "Delta E", "Time (s)")
        println("-"^78)
    end

    for iter in 1:1000
        t0 = time()

        build_total_J_matrix!(ws, orbitals)

        assemble_K_matrix!(ws, ws.K_mats[0], 0, orbitals)
        assemble_K_matrix!(ws, ws.K_mats[1], 1, orbitals)
        assemble_K_matrix!(ws, ws.K_mats[2], 2, orbitals)
        
        # --- Build Fock Matrices ---
        ws.F_s .= H_core_s .+ ws.J .- ws.K_mats[0]
        ws.F_p .= H_core_p .+ ws.J .- ws.K_mats[1]
        ws.F_d .= H_core_d .+ ws.J .- ws.K_mats[2]
        
        # --- Diagonalize Independent Blocks ---
        evals_fs, evecs_fs = eigen(Symmetric(ws.F_s[active_s, active_s]), ws.S[active_s, active_s])
        evals_fp, evecs_fp = eigen(Symmetric(ws.F_p[active_p, active_p]), ws.S[active_p, active_p])
        evals_fd, evecs_fd = eigen(Symmetric(ws.F_d[active_d, active_d]), ws.S[active_d, active_d])
        
        # --- Update Orbitals (With Mixing) ---
        
        # Update s-orbitals
        for (i, orb_idx) in enumerate(1:5)
            c_new = zeros(Float64, n); c_new[active_s] = evecs_fs[:, i]
            c_new ./= sqrt(dot(c_new, ws.S * c_new))
            orbitals[orb_idx].coeffs = MIXING * orbitals[orb_idx].coeffs + (1 - MIXING) * c_new
            orbitals[orb_idx].coeffs ./= sqrt(dot(orbitals[orb_idx].coeffs, ws.S * orbitals[orb_idx].coeffs))
            orbitals[orb_idx].energy = evals_fs[i]
        end

        # Update p-orbitals
        for (i, orb_idx) in enumerate(6:9)
            c_new = zeros(Float64, n); c_new[active_p] = evecs_fp[:, i]
            c_new ./= sqrt(dot(c_new, ws.S * c_new))
            orbitals[orb_idx].coeffs = MIXING * orbitals[orb_idx].coeffs + (1 - MIXING) * c_new
            orbitals[orb_idx].coeffs ./= sqrt(dot(orbitals[orb_idx].coeffs, ws.S * orbitals[orb_idx].coeffs))
            orbitals[orb_idx].energy = evals_fp[i]
        end

        # Update d-orbitals
        for (i, orb_idx) in enumerate(10:11)
            c_new = zeros(Float64, n); c_new[active_d] = evecs_fd[:, i]
            c_new ./= sqrt(dot(c_new, ws.S * c_new))
            orbitals[orb_idx].coeffs = MIXING * orbitals[orb_idx].coeffs + (1 - MIXING) * c_new
            orbitals[orb_idx].coeffs ./= sqrt(dot(orbitals[orb_idx].coeffs, ws.S * orbitals[orb_idx].coeffs))
            orbitals[orb_idx].energy = evals_fd[i]
        end

        # --- Compute Total Energy ---
        E_total = 0.0
        
        for orb in orbitals
            if orb.l == 0
                h_ii = dot(orb.coeffs, H_core_s * orb.coeffs)
            elseif orb.l == 1
                h_ii = dot(orb.coeffs, H_core_p * orb.coeffs)
            elseif orb.l == 2
                h_ii = dot(orb.coeffs, H_core_d * orb.coeffs)
            end
            E_total += (orb.occ / 2.0) * (h_ii + orb.energy)
        end
        
        delta = abs(E_total - E_old)
        elapsed = time() - t0
        
        if verbose
            @printf("%-4d | %14.8f | %10.2e | %8.4f\n", 
                    iter, E_total, delta, elapsed)
        end
        
        if delta < 1e-9
            if verbose 
                println("-"^78)
                println("Converged in $iter iterations.")
            end
            time2 = time()
            time_total = time2-time1
            println("Radio de confinamiento: $R_max a.u.")
            @printf("Energía final: %.6f Ha\n", E_total)
            println("Límite Hartree-Fock: -7232.1384 Ha")
            @printf("Tiempo elapsado: %.4f s\n", time_total)
            println("===== END =====")
            # --- Save the results ---
            filename = "xenon_results_R$(R_max).jld2"
            @save filename orbitals E_total R_max
            println("Saved orbitals and energy to $filename")
            
            # enforce_positive_phase!(orbitals)
            # plot_xenon_fischer_grid(ws.basis, orbitals, R_max)
            # plot_xenon_density(ws.basis, orbitals, R_max)
            break
        end
        
        E_old = E_total
    end
end
function plot_all_xenon_orbitals(basis::BSplineBasis, orbitals::Vector{Orbital}, R_max::Float64)
    # 1. Define a dense radial grid
    # Xenon's 5s and 5p extend quite far, so 20.0 a.u. is a good viewing window
    plot_max = min(R_max, 8.0) 
    r_points = collect(range(1e-5, plot_max, length=1000))
    
    # 2. Match the exact order of your Xenon orbitals array
    labels = ["1s", "2s", "3s", "4s", "5s", "2p", "3p", "4p", "5p", "3d", "4d"]
    
    # 3. Initialize the single main plot
    # Moving the legend to :outertopright prevents it from covering the wavefunctions
    p = plot(title="Radial Wavefunctions of Xenon (Z=54)", 
             xlabel="r (a.u.)", 
             ylabel="P(r)", 
             xlims=(0, plot_max),
             legend=:outertopright, 
             size=(850, 600),
             margin=5Plots.mm)
             
    # 4. Extract and plot each orbital
    for (i, orb) in enumerate(orbitals)
        psi_vals = evaluate_orbital(basis, orb.coeffs, r_points)
        
        # Use different line styles based on angular momentum to help tell them apart
        ls = if orb.l == 0
            :solid   # s-orbitals
        elseif orb.l == 1
            :dash    # p-orbitals
        else
            :dot     # d-orbitals
        end
        
        plot!(p, r_points, psi_vals, label=labels[i], lw=2, linestyle=ls)
    end
    
    # Display the final combined plot
    display(p)
    
    # Optional: Save it
    savefig(p, "xenon_all_orbitals.pdf")
end

function plot_xenon_fischer_grid(basis::BSplineBasis, orbitals::Vector{Orbital}, R_max::Float64)
    # 1. Recreate the Fischer grid limits
    # The textbook stops around r = 8.0 to 9.0, so we limit our viewing window.
    r_ticks = [0.0, 0.125, 0.5, 1.125, 2.0, 3.125, 4.5, 6.125, 8.0]
    sqrt_r_ticks = sqrt.(r_ticks)
    
    # Format labels to match the book (e.g., "1.12" instead of "1.125")
    tick_labels = ["0.0", "0.12", "0.5", "1.12", "2.0", "3.12", "4.5", "6.12", "8.0"]
    
    # 2. Define a very dense fine grid uniformly spaced in sqrt(r)
    # This ensures smooth curves when we plot against sqrt(r)
    sqrt_r_points = collect(range(0.0, sqrt(8.5), length=1500))
    r_points = sqrt_r_points .^ 2  # The actual physical r values to evaluate
    
    # Match the exact order of your Xenon orbitals array
    labels = ["1s", "2s", "3s", "4s", "5s", "2p", "3p", "4p", "5p", "3d", "4d"]
    
    # 3. Initialize the plot with the custom X-axis ticks
    p = plot(title="Hartree-Fock Radial Wavefunctions of Xenon", 
             xlabel="r (a.u.)", 
             ylabel="P(r)", 
             legend=:outertopright, 
             size=(850, 600),
             margin=5Plots.mm,
             xticks=(sqrt_r_ticks, tick_labels),
             xlims=(0, sqrt(8.5)))
             
    # 4. Extract and plot each orbital
    for (i, orb) in enumerate(orbitals)
        # Evaluate using the true 'r_points'
        psi_vals = evaluate_orbital(basis, orb.coeffs, r_points)
        
        # Use different line styles based on angular momentum
        ls = if orb.l == 0
            :solid
        elseif orb.l == 1
            :dash
        else
            :dot
        end
        
        # The key trick: we plot against `sqrt_r_points` instead of `r_points`!
        plot!(p, sqrt_r_points, psi_vals, label=labels[i], lw=2, linestyle=ls)
    end
    
    display(p)
    savefig(p, "xenon_all_orbitals.pdf")
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

function plot_xenon_density(basis::BSplineBasis, orbitals::Vector{Orbital}, R_max::Float64)
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
        
        # Add to total density: occupation * |P_nl(r)|^2
        for j in eachindex(total_density)
            total_density[j] += orb.occ * (psi_vals[j]^2)
        end
    end
    
    # 4. Create the plot
    p = plot(title="Total Radial Electron Density of Xenon", 
             xlabel="r (a.u.)", 
             ylabel="D(r)", 
             legend=false, 
             size=(850, 500),
             margin=5Plots.mm,
             xticks=(sqrt_r_ticks, tick_labels),
             xlims=(0, sqrt(8.5)))
             
    # Plot with a solid line and a subtle fill under the curve
    plot!(p, sqrt_r_points, total_density, lw=2, color=:black, fill=(0, 0.2, :blue))
    
    display(p)
    
    savefig(p, "xenon_total_density.pdf")
end

solve_xenon(20.0)
