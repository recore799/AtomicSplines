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


