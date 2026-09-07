using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))

using AtomicSplines
using LinearAlgebra
using Printf
using JLD2

# Import the modified solve_germanium_rohf function
include("germanium_rohf_vpol.jl")
include("np2_toy_ci.jl")
include("germanium_fine_structure.jl")

function run_calibration(alpha_range; r_c = 1.0)
    best_alpha = 0.0
    best_error = 1e9
    target_3P1 = 557.1
    target_3P2 = 1410.0
    
    println("Starting V_pol calibration for Germanium (Z=32.0)")
    println("Target splitting from 3P0:")
    println("  3P1 = $target_3P1 cm^-1")
    println("  3P2 = $target_3P2 cm^-1")
    println("-"^60)
    
    for alpha_d in alpha_range
        @printf("Testing alpha_d = %.3f (r_c = %.2f)...\n", alpha_d, r_c)
        
        # Run ROHF with vpol
        solve_germanium_rohf(30.0, alpha_d, r_c, "av"; verbose=false)
        
        # Run CI
        filename = @sprintf("germanium_rohf_results_av_R30.0_ad%.3f.jld2", alpha_d)
        
        # Load the newly saved effective potential and orbital to compute zeta
        data = load(filename)
        V_eff = data["V_eff"]
        P_4p = data["P_4p"]
        dense_grid = data["R_grid"]
        zeta_np = compute_zeta(dense_grid, V_eff, P_4p)
        println("Computed zeta_np = $zeta_np Hartrees")
        
        println("Running CI for alpha_d = $alpha_d")
        run_toy_ci("Germanium", filename, zeta_np) 
        println("-"^60)
    end
end

run_calibration([0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
