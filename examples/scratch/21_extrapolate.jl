using DelimitedFiles
using LinearAlgebra
using Printf

function calculate_cbs_limit(filename="ci_convergence.csv")
    if !isfile(filename)
        println("No se encontró el archivo $filename.")
        return
    end

    # Read data, skipping the header row
    data = readdlm(filename, ',', skipstart=1)
    
    L_max_vals = Float64.(data[:, 1])
    E_corr_vals = Float64.(data[:, 3])
    
    println("--- Datos de CI ---")
    for i in 1:length(L_max_vals)
        @printf("L_max = %d  |  E_corr = %.8f Ha\n", Int(L_max_vals[i]), E_corr_vals[i])
    end
    
    # Schwartz variable: X = 1 / (L + 1/2)^3
    X = 1.0 ./ (L_max_vals .+ 0.5).^3
    
    # Least squares fit: E_corr = E_inf + A * X
    # M * [E_inf; A] = E_corr_vals
    M = hcat(ones(length(X)), X)
    coeffs = M \ E_corr_vals
    
    E_inf = coeffs[1]
    A_coeff = coeffs[2]
    
    println("\n--- Límite Extrapolado (CBS) ---")
    @printf("E_corr (L -> ∞) : %.6f Ha\n", E_inf)
    @printf("Coeficiente A   : %.6f\n", A_coeff)
    println("--------------------------------")
end

calculate_cbs_limit()
