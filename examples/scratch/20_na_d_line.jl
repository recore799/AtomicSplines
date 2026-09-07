using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using Plots
using JLD2
using LinearAlgebra
using Printf

function compute_sodium_d_line(data_file::String)
    # 1. Load the data
    data = load(data_file)
    
    R_mat = data["R_mat"]
    evals_fs = data["evals_fs"]
    evecs_fs = data["evecs_fs"]
    evals_fp = data["evals_fp"]
    evecs_fp = data["evecs_fp"]
    active_s = data["active_s"]
    active_p = data["active_p"]
    n = data["num_splines"]

    # 2. Extract the 3s Valence Ground State (Index 3 of the s-block)
    E_3s = evals_fs[3]
    c_3s = zeros(Float64, n)
    c_3s[active_s] .= evecs_fs[:, 3]

    # 3. Extract the 3p Excited State (Index 2 of the p-block)
    E_3p = evals_fp[2]
    c_3p = zeros(Float64, n)
    c_3p[active_p] .= evecs_fp[:, 2]

    # 4. Calculate Transition Energy
    dE_Ha = E_3p - E_3s
    dE_eV = dE_Ha * 27.211386
    wavelength_nm = 1239.84193 / dE_eV

    # 5. Calculate the Radial Dipole Matrix Element <3s | r | 3p>
    radial_dipole = dot(c_3s, R_mat * c_3p)
    dipole_sq = radial_dipole^2

    # 6. Print the Results
    println("\n=== Sodium D-Line (3p -> 3s) ===")
    @printf("Energy 3s : %12.6f Ha\n", E_3s)
    @printf("Energy 3p : %12.6f Ha\n", E_3p)
    println("-"^34)
    @printf("Delta E   : %12.6f Ha\n", dE_Ha)
    @printf("Delta E   : %12.6f eV\n", dE_eV)
    @printf("Wavelength: %12.2f nm\n", wavelength_nm)
    @printf("Dipole Sq : %12.6e\n", dipole_sq)

    # 1s is index 1, 2s is index 2 in the s-block
    c_1s = zeros(Float64, n)
    c_1s[active_s] .= evecs_fs[:, 1]
    
    c_2s = zeros(Float64, n)
    c_2s[active_s] .= evecs_fs[:, 2]
    
    # 2p is index 1 in the p-block
    c_2p = zeros(Float64, n)
    c_2p[active_p] .= evecs_fp[:, 1]

    # Recreate the exact basis used for the SCF run
    R_max = 50.0 
    N_elems = 300
    basis = generate_basis(R_max, N_elems, Val(7), γ=2.5)

    # Generate the physically complete figure!
    plot_total_sodium_penetration(basis, c_1s, c_2s, c_2p, c_3s, c_3p)

end
function plot_total_sodium_penetration(basis::BSplineBasis, c_1s, c_2s, c_2p, c_3s, c_3p)
    # 1. Grid setup (using the sqrt scale)
    sqrt_r_points = collect(range(0.0, sqrt(15.0), length=1500))
    r_points = sqrt_r_points .^ 2  
    
    # 2. Evaluate all wavefunctions
    psi_1s = evaluate_orbital(basis, c_1s, r_points)
    psi_2s = evaluate_orbital(basis, c_2s, r_points)
    psi_2p = evaluate_orbital(basis, c_2p, r_points)
    psi_3s = evaluate_orbital(basis, c_3s, r_points)
    psi_3p = evaluate_orbital(basis, c_3p, r_points)
    
    # 3. Calculate Total Core Density
    # Na+ core: 1s², 2s², 2p⁶
    core_dens = 2.0 .* (psi_1s .^ 2) .+ 
                2.0 .* (psi_2s .^ 2) .+ 
                6.0 .* (psi_2p .^ 2)

    dens_3s = psi_3s .^ 2
    dens_3p = psi_3p .^ 2

    # 4. Create the plot
    r_ticks = [0.0, 0.5, 2.0, 4.5, 8.0, 12.5]
    sqrt_r_ticks = sqrt.(r_ticks)
    tick_labels = ["0.0", "0.5", "2.0", "4.5", "8.0", "12.5"]
    
    p = plot(title="Sodium: Valence Penetration into Total Core", 
             xlabel="r (a.u.) [sqrt scale]", 
             ylabel="Radial Density P²(r)", 
             legend=:topright, 
             size=(850, 600),
             margin=5*Plots.mm,
             xticks=(sqrt_r_ticks, tick_labels),
             xlims=(0, sqrt(15.0)))
             
    # 5. Plot the Total Core (Shaded)
    plot!(p, sqrt_r_points, core_dens, 
          label="Total Core (1s² + 2s² + 2p⁶)", 
          fill=true, fillalpha=0.3, color=:gray, 
          linewidth=1)

    # 6. Plot the Valence States (Magnified for visibility)
    # mag_factor = 200.0
    mag_factor = 20.0
    
    plot!(p, sqrt_r_points, dens_3s .* mag_factor, 
          label="3s Valence (x20)", color=:blue, lw=1, linestyle=:solid)
          
    plot!(p, sqrt_r_points, dens_3p .* mag_factor, 
          label="3p Valence (x20)", color=:red, lw=1, linestyle=:solid)
    
    display(p)
    savefig(p, "sodium_penetration.pdf")
end

compute_sodium_d_line("sodium_core_results.jld2")
