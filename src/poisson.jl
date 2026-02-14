"""
    solve_poisson_tensor(ws::SolverWorkspace{K}, orbital_coeffs)

Solves Poisson using pre-computed interaction tensors.
Replaces the quadrature loop with a tensor contraction: b_k = c^T * W_k * c
"""
function solve_poisson_tensor(ws::SolverWorkspace{K}, orbital_coeffs::Vector{Float64}) where {K}
    basis = ws.basis
    n = basis.num_splines
    
    # 1. Initialize RHS vector b
    b = zeros(Float64, n)
    
    # 2. Tensor Contraction Loop (Iterate over Elements)
    for i in 1:(length(basis.knots) - 1)
        # Skip if tensor wasn't built (e.g., empty interval)
        if !isassigned(ws.interaction_tensors, i); continue; end
        
        # Get the pre-computed physics for this element
        W_local = ws.interaction_tensors[i] 
        first_global = i - K + 1

        # A. Extract Local Coefficients (c_a)
        # We need the density c_a * c_b in this element
        c_local = zeros(Float64, K)
        for idx in 1:K
            g_idx = first_global + idx - 1
            if g_idx >= 1 && g_idx <= n
                c_local[idx] = orbital_coeffs[g_idx]
            end
        end

        # B. Contract: b[k] += sum_{a,b} ( c[a] * c[b] * W[k,a,b] )
        # This calculates the projection of the density onto the basis B_k
        for k in 1:K
            g_k = first_global + k - 1
            if g_k < 1 || g_k > n; continue; end
            
            val_k = 0.0
            for a in 1:K
                c_a = c_local[a]
                if abs(c_a) < 1e-15; continue; end
                
                # W is symmetric in a,b. We can optimize the inner loop.
                # Term (a,a)
                val_k += c_a * c_a * W_local[k, a, a]
                
                # Terms (a,b) and (b,a) are equal
                for b in (a+1):K
                    val_k += 2.0 * c_a * c_local[b] * W_local[k, a, b]
                end
            end
            b[g_k] += val_k
        end
    end

    # 3. Apply Boundary Conditions (Same as before)
    active = 2:(n-1)
    y_boundary = 1.0 
    
    # Shift RHS for Dirichlet BCs
    boundary_shift = ws.T[active, n] 
    boundary_shift .*= (2.0 * y_boundary) # K_stiff = 2*T
    
    b_inner = b[active] 
    b_inner .-= boundary_shift
    
    # 4. Fast Solve
    y_inner = ws.poisson_fact \ b_inner
    
    y_full = zeros(Float64, n)
    y_full[active] = y_inner
    y_full[n] = y_boundary
    
    return y_full
end
"""
    solve_poisson_fast(ws::SolverWorkspace{K}, orbital_coeffs)

Solves the radial Poisson equation -y'' = u^2/r for the Hartree potential.
Uses the pre-computed Cholesky factorization from the workspace.
"""
function solve_poisson_fast(ws::SolverWorkspace{K}, orbital_coeffs::Vector{Float64}) where {K}
    basis = ws.basis
    n = basis.num_splines
    
    # 1. Reset RHS vector b (re-use allocated memory if possible, but creating a new one is cheap O(N))
    b = zeros(Float64, n)
    
    # 2. Assembly Loop (Allocation Free)
    #    Iterate over elements (intervals)
    for i in 1:(length(basis.knots) - 1)
        t_a = basis.knots[i]; t_b = basis.knots[i+1]
        if t_a == t_b; continue; end
        
        mid = (t_a + t_b)/2; scale = (t_b - t_a)/2
        first_global = i - K + 1

        # Iterate Quadrature Points (Pre-cached in basis)
        for q in 1:length(basis.gl_nodes)
            r = scale * basis.gl_nodes[q] + mid
            w = scale * basis.gl_weights[q]

            # A. Kernel Call (Unrolled for K)
            eval_bspline_kernel!(ws.scratch_vals, ws.scratch_derivs, Val(true), Val(false), 
                                 i, r, basis.knots, Val(K))

            # B. Reconstruct Radial Function u(r) locally
            u_val = 0.0
            for idx in 1:K
                g_idx = first_global + idx - 1
                if g_idx >= 1 && g_idx <= n
                    u_val += orbital_coeffs[g_idx] * ws.scratch_vals[idx]
                end
            end
            
            # C. Source Term: u^2 / r
            #    Handle singularity safely
            if r > 1e-10
                source = (u_val^2) / r
            else
                source = 0.0
            end

            # D. Scatter into RHS vector b
            for a in 1:K
                g_a = first_global + a - 1
                if g_a >= 1 && g_a <= n
                    b[g_a] += w * ws.scratch_vals[a] * source
                end
            end
        end
    end

    # 3. Apply Boundary Conditions & Solve
    #    BCs: y(0)=0 (implicit via slice), y(R) = 1.0 (Total Charge)
    
    # Define active region (internal nodes)
    active = 2:(n-1)
    y_boundary = 1.0 
    
    # Retrieve Stiffness Matrix Column 'n' for the boundary shift
    # K_stiff = 2 * T. We can compute the column on the fly or grab it from T.
    # Shift = K[inner, n] * y_b = 2 * T[inner, n] * 1.0
    
    # Note: accessing column 'n' of a dense matrix is fast.
    boundary_shift = ws.T[active, n] 
    boundary_shift .*= (2.0 * y_boundary) # Scale by 2 (Kinetic -> Stiffness)
    
    # Adjust RHS
    b_inner = b[active] 
    b_inner .-= boundary_shift
    
    # 4. Fast Solve (Cholesky)
    #    ws.poisson_fact is the pre-computed factorization of 2*T[active, active]
    y_inner = ws.poisson_fact \ b_inner
    
    # 5. Reconstruct Full Vector
    y_full = zeros(Float64, n)
    y_full[active] = y_inner
    y_full[n] = y_boundary
    
    return y_full
end


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
            if r > 1e-1
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




