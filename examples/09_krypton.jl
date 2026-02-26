using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using LinearAlgebra
using Plots
using Printf
using JLD2

function solve_krypton(R_max; verbose::Bool=true)
    println("=== Krypton (Z=36) ===")

    time1 = time()
    
    N_elems = 100
    Z = 36.0
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
        Orbital(2, 1, 6.0),  # 2p (idx 5)
        Orbital(3, 1, 6.0),  # 3p (idx 6)
        Orbital(4, 1, 6.0),  # 4p (idx 7)
        Orbital(3, 2, 10.0)  # 3d (idx 8)
    ]
    
    # --- Core Hamiltonians ---
    H_core_s = ws.T + ws.V
    H_core_p = ws.T + ws.V + ws.V2 
    H_core_d = ws.T + ws.V + 3.0 * ws.V2  # l(l+1)/2 = 6/2 = 3
    
    # --- Initial Guess ---
    evals_s, evecs_s = eigen(Symmetric(H_core_s[active_s, active_s]), ws.S[active_s, active_s])
    evals_p, evecs_p = eigen(Symmetric(H_core_p[active_p, active_p]), ws.S[active_p, active_p])
    evals_d, evecs_d = eigen(Symmetric(H_core_d[active_d, active_d]), ws.S[active_d, active_d])
    
    # Assign initial guesses based on principal quantum number sorting
    for (i, orb_idx) in enumerate(1:4) # 1s, 2s, 3s, 4s
        orbitals[orb_idx].coeffs = zeros(Float64, n)
        orbitals[orb_idx].coeffs[active_s] = evecs_s[:, i]
    end
    for (i, orb_idx) in enumerate(5:7) # 2p, 3p, 4p
        orbitals[orb_idx].coeffs = zeros(Float64, n)
        orbitals[orb_idx].coeffs[active_p] = evecs_p[:, i]
    end
    # 3d
    orbitals[8].coeffs = zeros(Float64, n)
    orbitals[8].coeffs[active_d] = evecs_d[:, 1]
    
    # Normalize initial guesses
    for orb in orbitals
        orb.coeffs ./= sqrt(dot(orb.coeffs, ws.S * orb.coeffs))
    end

    # --- SCF Loop ---
    E_old = 0.0
    MIXING = 0.3 # Krypton often needs higher damping to prevent oscillations
    
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
        for (i, orb_idx) in enumerate(1:4)
            c_new = zeros(Float64, n); c_new[active_s] = evecs_fs[:, i]
            c_new ./= sqrt(dot(c_new, ws.S * c_new))
            orbitals[orb_idx].coeffs = MIXING * orbitals[orb_idx].coeffs + (1 - MIXING) * c_new
            orbitals[orb_idx].coeffs ./= sqrt(dot(orbitals[orb_idx].coeffs, ws.S * orbitals[orb_idx].coeffs))
            orbitals[orb_idx].energy = evals_fs[i]
        end

        # Update p-orbitals
        for (i, orb_idx) in enumerate(5:7)
            c_new = zeros(Float64, n); c_new[active_p] = evecs_fp[:, i]
            c_new ./= sqrt(dot(c_new, ws.S * c_new))
            orbitals[orb_idx].coeffs = MIXING * orbitals[orb_idx].coeffs + (1 - MIXING) * c_new
            orbitals[orb_idx].coeffs ./= sqrt(dot(orbitals[orb_idx].coeffs, ws.S * orbitals[orb_idx].coeffs))
            orbitals[orb_idx].energy = evals_fp[i]
        end

        # Update d-orbital
        c_new_d = zeros(Float64, n); c_new_d[active_d] = evecs_fd[:, 1]
        c_new_d ./= sqrt(dot(c_new_d, ws.S * c_new_d))
        orbitals[8].coeffs = MIXING * orbitals[8].coeffs + (1 - MIXING) * c_new_d
        orbitals[8].coeffs ./= sqrt(dot(orbitals[8].coeffs, ws.S * orbitals[8].coeffs))
        orbitals[8].energy = evals_fd[1]

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
            println("Límite Hartree-Fock: -2752.0550 Ha")
            @printf("Tiempo elapsado: %.4f s\n", time_total)
            println("===== END =====")

            # --- Save the results ---
            filename = "krypton_results_R$(R_max).jld2"
            @save filename orbitals E_total R_max
            println("Saved orbitals and energy to $filename")
            
            # plot_krypton_orbitals(ws.basis, orbitals, R_max)
            break
        end
        
        E_old = E_total
    end
end

function plot_krypton_orbitals(basis::BSplineBasis, orbitals::Vector{Orbital}, R_max::Float64)
    # 1. Define a dense radial grid for smooth plotting
    # We focus on the inner region (e.g., up to 10 a.u.) where most core density lives,
    # but you can extend it to R_max if you want to see the long tails.
    plot_max = min(R_max, 8.0) 
    r_points = collect(range(1e-5, plot_max, length=1000))
    
    # 2. Extract evaluated functions
    # Assuming your Krypton orbitals array is ordered: 
    # [1s, 2s, 3s, 4s, 2p, 3p, 4p, 3d]
    labels = ["1s", "2s", "3s", "4s", "2p", "3p", "4p", "3d"]
    psi_dict = Dict{String, Vector{Float64}}()
    
    for (i, orb) in enumerate(orbitals)
        # Your evaluate_orbital function handles the B-spline contraction
        psi_dict[labels[i]] = evaluate_orbital(basis, orb.coeffs, r_points)
    end

    # 3. Create subplots grouped by angular momentum
    
    # --- Plot s-orbitals (l=0) ---
    p_s = plot(title="s-orbitals (l=0)", xlabel="r (a.u.)", ylabel="P(r)", xlims=(0, plot_max))
    plot!(p_s, r_points, psi_dict["1s"], label="1s", lw=2)
    plot!(p_s, r_points, psi_dict["2s"], label="2s", lw=2)
    plot!(p_s, r_points, psi_dict["3s"], label="3s", lw=2)
    plot!(p_s, r_points, psi_dict["4s"], label="4s", lw=2)
    
    # --- Plot p-orbitals (l=1) ---
    p_p = plot(title="p-orbitals (l=1)", xlabel="r (a.u.)", ylabel="P(r)", xlims=(0, plot_max))
    plot!(p_p, r_points, psi_dict["2p"], label="2p", lw=2)
    plot!(p_p, r_points, psi_dict["3p"], label="3p", lw=2)
    plot!(p_p, r_points, psi_dict["4p"], label="4p", lw=2)
    
    # --- Plot d-orbitals (l=2) ---
    p_d = plot(title="d-orbitals (l=2)", xlabel="r (a.u.)", ylabel="P(r)", xlims=(0, plot_max))
    plot!(p_d, r_points, psi_dict["3d"], label="3d", lw=2, color=:purple)
    
    # 4. Combine into a single layout
    final_plot = plot(p_s, p_p, p_d, layout=(3, 1), size=(800, 900), margin=5Plots.mm)
    
    # Display the plot
    display(final_plot)
    
    # Optional: Save it to a file
    # savefig(final_plot, "krypton_radial_functions.png")
end

solve_krypton(20.0)
