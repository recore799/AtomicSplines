function init_scf_workspace(basis::BSplineBasis{K}, Z::Float64) where {K}
    n = basis.num_splines
    
    # Allocate Matrices
    T = zeros(Float64, n, n)
    V = zeros(Float64, n, n)
    S = zeros(Float64, n, n)
    
    # Scratch buffers (Exact size K)
    scratch_v = zeros(Float64, K)
    scratch_d = zeros(Float64, K)
    
    for i in 1:(length(basis.knots) - 1)
        t_a = basis.knots[i]; t_b = basis.knots[i+1]
        if t_a == t_b; continue; end
        
        mid = (t_a + t_b)/2; scale = (t_b - t_a)/2
        first_global = i - K + 1

        for q in 1:length(basis.gl_nodes)
            r = scale * basis.gl_nodes[q] + mid
            w = scale * basis.gl_weights[q]

            # Pass Val{K} to kernel
            eval_bspline_kernel!(scratch_v, scratch_d, Val(true), Val(true), 
                                 i, r, basis.knots, Val(K))

            # Block Scatter
            for a in 1:K
                g_a = first_global + a - 1
                if g_a < 1 || g_a > n; continue; end
                
                Na = scratch_v[a]; dNa = scratch_d[a]
                
                for b in a:K # Symmetry optimization
                    g_b = first_global + b - 1
                    if g_b < 1 || g_b > n; continue; end
                    
                    Nb = scratch_v[b]; dNb = scratch_d[b]
                    
                    S[g_a, g_b] += w * Na * Nb
                    T[g_a, g_b] += w * 0.5 * dNa * dNb
                    if r > 1e-12
                        V[g_a, g_b] += w * (-Z / r) * Na * Nb
                    end
                end
            end
        end
    end
    
    # Symmetrize
    S = Symmetric(S); T = Symmetric(T); V = Symmetric(V)

    # Build the Interaction Tensors
    tensors = build_interaction_tensors(basis)

    # --- PRE-COMPUTE POISSON (k=0) ---
    # Stiffness K = 2*T. Factorize inner block.
    K_mat = 2.0 .* T
    active = 2:(n-1)
    
    # This is the "Cholesky Trick" you wanted
    # We store the factorization, not the matrix!
    poisson_fact = cholesky(K_mat[active, active])
    
    return SolverWorkspace{K}(basis, Matrix(T), Matrix(V), Matrix(S), 
                              tensors, poisson_fact, scratch_v, scratch_d)
end


"""
Builds T (Kinetic), V_nuc (Nuclear), and S (Overlap) in a single pass.
Iterates over ELEMENTS (intervals), not functions.
"""
function assemble_core(basis, Z)
    n = basis.num_splines
    k = basis.order
    
    # Use dense matrices for N < 500. 
    T = zeros(Float64, n, n)
    V = zeros(Float64, n, n)
    S = zeros(Float64, n, n)
    
    scratch = BSplineScratch(Vector{Float64}(undef, k), Vector{Float64}(undef, k))
    
    # Gauss-Legendre Quadrature
    gl_p, gl_w = gausslegendre(k + 2) 

    # Iterate over knot intervals (elements)
    # The active range of splines is roughly [1, n+k] in the knot vector
    for i in 1:(length(basis.knots) - 1)
        t_a = basis.knots[i]; t_b = basis.knots[i+1]
        if t_a == t_b; continue; end

        # Map bounds
        mid = (t_a + t_b) / 2
        scale = (t_b - t_a) / 2
        
        # Identify active B-splines: B_{i-k+1} ... B_{i} are non-zero here
        first_global = i - k + 1

        for q in 1:length(gl_p)
            r = scale * gl_p[q] + mid
            w = scale * gl_w[q]

            # Kernel Call: Get Vals AND Derivs
            eval_bspline_kernel!(scratch.vals, scratch.derivs, Val(true), Val(true), 
                                 i, k, r, basis.knots)

            # Add to matrix blocks
            for a in 1:k
                g_a = first_global + a - 1
                if g_a < 1 || g_a > n; continue; end
                
                # Pre-fetch for inner loop
                Na  = scratch.vals[a]
                dNa = scratch.derivs[a]

                for b in 1:k
                    g_b = first_global + b - 1
                    if g_b < 1 || g_b > n; continue; end
                    
                    Nb  = scratch.vals[b]
                    dNb = scratch.derivs[b]
                    
                    # Overlap: <a|b>
                    S[g_a, g_b] += w * Na * Nb
                    
                    # Kinetic: 0.5 * <da|db>
                    T[g_a, g_b] += w * 0.5 * dNa * dNb
                    
                    # Nuclear: -Z * <a|1/r|b>
                    # Avoid singularity at r=0 (though GL points shouldn't hit 0)
                    if r > 1e-12
                        V[g_a, g_b] += w * (-Z / r) * Na * Nb
                    end
                end
            end
        end
    end
    
    return Symmetric(T), Symmetric(V), Symmetric(S)
end

function assemble_J_matrix_tensor(ws::SolverWorkspace{K}, y_coeffs) where {K}
    n = ws.basis.num_splines
    J_mat = zeros(Float64, n, n)
    
    for i in 1:(length(ws.basis.knots) - 1)
        # Retrieve pre-computed physics
        if !isassigned(ws.interaction_tensors, i); continue; end
        W_local = ws.interaction_tensors[i]
        
        first_global = i - K + 1
        
        # 1. Get local y coefficients
        # "Reconstruct" is now just grabbing indices
        y_local = zeros(Float64, K)
        for idx in 1:K
            g_idx = first_global + idx - 1
            if g_idx >= 1 && g_idx <= n
                y_local[idx] = y_coeffs[g_idx]
            end
        end
        
        # 2. Contract Tensor: J_ab = sum_k ( y_k * W_kab )
        # This is effectively adding weighted matrices
        for k in 1:K
            weight = y_local[k]
            if abs(weight) < 1e-12; continue; end
            
            # Add the slice W[k, :, :] scaled by y[k]
            for a in 1:K
                g_a = first_global + a - 1
                if g_a < 1 || g_a > n; continue; end
                
                for b in 1:K
                    g_b = first_global + b - 1
                    if g_b < 1 || g_b > n; continue; end
                    
                    J_mat[g_a, g_b] += weight * W_local[k, a, b]
                end
            end
        end
    end
    
    return Symmetric(J_mat)
end

function assemble_J_matrix_param(ws::SolverWorkspace{K}, y_coeffs) where {K}
    basis = ws.basis
    n = basis.num_splines

    J_mat = zeros(Float64, n, n)

    for i in 1:(length(basis.knots) - 1)
        t_a = basis.knots[i]; t_b = basis.knots[i+1]
        if t_a == t_b; continue; end

        mid = (t_a + t_b)/2; scale = (t_b - t_a)/2
        first_global = i - K + 1

        for q in 1:length(basis.gl_nodes)
            r = scale * basis.gl_nodes[q] + mid
            w = scale * basis.gl_weights[q]

            eval_bspline_kernel!(ws.scratch_vals, ws.scratch_derivs, Val(true), Val(true),
                                 i, r, basis.knots, Val(K))
        
            # Compute y(r) = r * V_H(r)
            # Actually, usually y_coeffs ARE the density projection coeffs.
            # Let's assume y_coeffs represents the solution to Poisson: U(r) = y(r)/r
            # So here we reconstruct y(r)
            y_val = 0.0
            for idx in 1:K
                g_idx = first_global + idx - 1
                if g_idx >= 1 && g_idx <= n
                    y_val += y_coeffs[g_idx] * ws.scratch_vals[idx]
                end
            end
            
            # Potential V_H = y(r) / r
            if r < 1e-12; continue; end # Avoid singularity
            V_H = y_val / r
            
            # Fill Matrix J_ab = <a | V_H | b>
            for a in 1:K
                g_a = first_global + a - 1
                if g_a < 1 || g_a > n; continue; end
                
                term_a = ws.scratch_vals[a] * w * V_H
                
                for b in 1:K
                    g_b = first_global + b - 1
                    if g_b < 1 || g_b > n; continue; end
                
                    J_mat[g_a, g_b] += term_a * ws.scratch_vals[b]
                end
            end
        end
    end
    return Symmetric(J_mat)
end

"""
Builds the Hartree Potential Matrix J.
Calculates density y(r) on the fly without extra searches.
"""
function assemble_J_matrix(basis, y_coeffs)
    n = basis.num_splines; k = basis.order
    J_mat = zeros(Float64, n, n)
    
    scratch = BSplineScratch(Vector{Float64}(undef, k), Vector{Float64}(undef, k))
    gl_p, gl_w = gausslegendre(k + 2)

    for i in 1:(length(basis.knots) - 1)
        t_a = basis.knots[i]; t_b = basis.knots[i+1]
        if t_a == t_b; continue; end
        
        first_global = i - k + 1
        mid = (t_a + t_b)/2; scale = (t_b - t_a)/2

        for q in 1:length(gl_p)
            r = scale * gl_p[q] + mid
            w = scale * gl_w[q]

            # Need Values only
            eval_bspline_kernel!(scratch.vals, scratch.derivs, Val(true), Val(false), 
                                 i, k, r, basis.knots)

            # Compute y(r) = r * V_H(r)
            # Actually, usually y_coeffs ARE the density projection coeffs.
            # Let's assume y_coeffs represents the solution to Poisson: U(r) = y(r)/r
            # So here we reconstruct y(r)
            y_val = 0.0
            for idx in 1:k
                g_idx = first_global + idx - 1
                if g_idx >= 1 && g_idx <= n
                    y_val += y_coeffs[g_idx] * scratch.vals[idx]
                end
            end
            
            # Potential V_H = y(r) / r
            if r < 1e-12; continue; end # Avoid singularity
            V_H = y_val / r
            
            # Fill Matrix J_ab = <a | V_H | b>
            for a in 1:k
                g_a = first_global + a - 1
                if g_a < 1 || g_a > n; continue; end
                
                term_a = scratch.vals[a] * w * V_H
                
                for b in 1:k
                    g_b = first_global + b - 1
                    if g_b < 1 || g_b > n; continue; end
                    
                    J_mat[g_a, g_b] += term_a * scratch.vals[b]
                end
            end
        end
    end
    return Symmetric(J_mat)
end


function build_interaction_tensors(basis::BSplineBasis{K}) where {K}
    tensors = Vector{Array{Float64, 3}}(undef, length(basis.knots) - 1)

    # Scratch for kernel
    vals = zeros(Float64, K); derivs = zeros(Float64, K)

    for i in 1:(length(basis.knots) - 1)
        t_a = basis.knots[i]; t_b = basis.knots[i+1]
        if t_a == t_b; continue; end
        mid = (t_a + t_b)/2; scale = (t_b - t_a)/2

        # Local tensor for this element
        W_local = zeros(Float64, K, K, K)

        for q in 1:length(basis.gl_nodes)
            r = scale * basis.gl_nodes[q] + mid
            w = scale * basis.gl_weights[q]

            # Singularity check
            inv_r = (r > 1e-10) ? 1.0/r : 0.0

            # Eval basis
            eval_bspline_kernel!(vals, derivs, Val(true), Val(false), i, r, basis.knots, Val(K))

            for k_idx in 1:K
                val_k = vals[k_idx] * w * inv_r
                for a in 1:K
                    val_ka = val_k * vals[a]
                    for b in a:K # Symmetric in a,b
                        W_local[k_idx, a, b] += val_ka * vals[b]
                        W_local[k_idx, b, a] = W_local[k_idx, a, b]
                    end
                end
            end
        end
        tensors[i] = W_local
    end
    return tensors
end
