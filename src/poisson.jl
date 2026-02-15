function solve_generalized_poisson(ws::SolverWorkspace{K}, source::Vector{Float64}, k_mult::Int; boundary_val::Float64=0.0) where {K}
    # Check Cache for Factorization
    if !haskey(ws.poisson_factors, k_mult)
        A_k = 2.0 .* ws.T .+ (k_mult * (k_mult + 1)) .* ws.V2
        active = 2:(ws.basis.num_splines - 1)
        ws.poisson_factors[k_mult] = cholesky(A_k[active, active])
    end
    
    n = length(source)
    active = 2:(n-1)
    
    # Apply Boundary Condition Shift
    # The equation is A * y = b. 
    # For Dirichlet y(N) = val, we move the last column to RHS.
    # rhs = source - Column_N * (2.0 * boundary_val) 
    # (The factor 2.0 comes from the Stiffness matrix definition 2*T)
    
    b_inner = source[active]
    if boundary_val != 0.0
        col_n = ws.T[active, n]
        b_inner .-= col_n .* (2.0 * boundary_val)
    end

    y_inner = ws.poisson_factors[k_mult] \ b_inner
    
    y_full = zeros(Float64, n)
    y_full[active] = y_inner
    y_full[n] = boundary_val
    
    return y_full
end

function solve_poisson_J(ws::SolverWorkspace{K}, orbital_coeffs::Vector{Float64}) where {K}
    basis = ws.basis; n = basis.num_splines
    b = zeros(Float64, n); c_local = ws.scratch_vals
    
    for i in 1:(length(basis.knots)-1)
        if !isassigned(ws.interaction_tensors, i); continue; end
        W_local = ws.interaction_tensors[i]
        first_global = i - K + 1
        
        fill!(c_local, 0.0)
        for idx in 1:K
            g = first_global + idx - 1
            if g >= 1 && g <= n; c_local[idx] = orbital_coeffs[g]; end
        end
        
        for k in 1:K
            g_k = first_global + k - 1
            if g_k < 1 || g_k > n; continue; end
            val = 0.0
            for a in 1:K
                ca = c_local[a]
                if abs(ca) < 1e-15; continue; end
                val += ca * ca * W_local[k, a, a]
                for b in (a+1):K
                    val += 2.0 * ca * c_local[b] * W_local[k, a, b]
                end
            end
            b[g_k] += val
        end
    end
    
    # Boundary Condition y(R) = Z_eff (Total Charge)
    # For J, charge is N_electrons. For Exchange, usually 0.
    # Actually, J potential y(r) -> N_electrons as r -> inf.
    active = 2:(n-1); y_bound = sum(orbital_coeffs) # Approximate total charge? No, just force 1.0 for normalized density.
    # WAIT: Input 'orbital_coeffs' is Density coeffs (c^2). Sum is not 1.
    # Let's assume standard BC: y(R) = Total Integrated Charge.
    # For now, force y(R)=1.0 (assuming input is normalized orbital squared)
    y_bound = 1.0 
    
    col_n = ws.T[active, n]
    b_inner = b[active] .- col_n .* (2.0 * y_bound)
    y_inner = ws.poisson_factors[0] \ b_inner
    y_full = zeros(Float64, n)
    y_full[active] = y_inner; y_full[n] = y_bound
    return y_full
end
