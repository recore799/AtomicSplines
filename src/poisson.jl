function solve_poisson_potential(basis, orbital_coeffs, T_matrix)
    n = basis.num_splines
    k = basis.order
    
    # The RHS vector b
    b = zeros(Float64, n)
    
    scratch = BSplineScratch(Vector{Float64}(undef, k), Vector{Float64}(undef, k))
    # Use slightly higher quadrature for the Coulomb term (degree ~3k)
    gl_p, gl_w = gausslegendre(k + 4)

    for i in 1:(length(basis.knots) - 1)
        t_a = basis.knots[i]; t_b = basis.knots[i+1]
        if t_a == t_b; continue; end
        
        first_global = i - k + 1
        mid = (t_a + t_b)/2; scale = (t_b - t_a)/2

        for q in 1:length(gl_p)
            r = scale * gl_p[q] + mid
            w = scale * gl_w[q]
            
            # 1. Eval basis
            eval_bspline_kernel!(scratch.vals, scratch.derivs, Val(true), Val(false), 
                                 i, k, r, basis.knots)
            
            # 2. Reconstruct radial function u(r)
            u_val = 0.0
            for idx in 1:k
                g_idx = first_global + idx - 1
                if g_idx >= 1 && g_idx <= n
                    u_val += orbital_coeffs[g_idx] * scratch.vals[idx]
                end
            end
            
            # 3. Source Term: -y'' = u^2 / r
            #    We solve K*y = b. K is positive definite (integral of B' B').
            #    The operator is -d^2/dr^2. So LHS is positive.
            #    RHS must be positive: u^2/r.
            if r > 1e-12
                source = u_val^2 / r
            else
                source = 0.0 
            end
            
            # 4. Fill b vector
            for a in 1:k
                g_a = first_global + a - 1
                if g_a < 1 || g_a > n; continue; end
                b[g_a] += w * scratch.vals[a] * source
            end
        end
    end
    
    # 5. Solve K * y = b
    #    K = integral(B' B') = 2 * T_matrix (stiffness)
    K = 2.0 .* T_matrix
    
    # --- BOUNDARY CONDITIONS ---
    # y(0) = 0 (Potential is finite, y=rV -> 0)
    # y(R) = 1 (Total Charge Enclosed = 1 electron)
    #
    # We solve for inner nodes [2:n-1].
    # The equation for node (n-1) involves y_n:
    # K[n-1, ...] * y_... + K[n-1, n] * y_n = b[n-1]
    # Move known boundary term to RHS:
    
    inner = 2:(n-1)
    y_boundary = 1.0 
    
    # Adjust RHS
    b_inner = b[inner] - K[inner, n] * y_boundary
    
    # Solve
    y_inner = K[inner, inner] \ b_inner
    
    # Reconstruct full vector
    y_full = zeros(Float64, n)
    y_full[inner] = y_inner
    y_full[n] = y_boundary # Enforce y(R) = 1
    
    return y_full
end




