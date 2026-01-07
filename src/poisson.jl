function solve_poisson_potential(basis, P_coeffs, T_matrix)
    n = basis.num_splines; k = basis.order
    
    F = zeros(n)
    gl_p, gl_w = gausslegendre(k + 2)
    
    for i in 1:n
        s_knot = i; e_knot = i+k
        for interval in s_knot:(e_knot-1)
            ta = basis.knots[interval]; tb = basis.knots[interval+1]; if ta==tb continue end
            mid = (ta+tb)/2; scale = (tb-ta)/2
            for q in 1:length(gl_p)
                r = scale*gl_p[q] + mid; w = scale*gl_w[q]
                P_val = eval_expansion(P_coeffs, basis, r)
                density_source = (P_val^2) / r
                F[i] += w * bspline(i, k, r, basis.knots) * density_source
            end
        end
    end
    
    K = 2.0 .* T_matrix
    
    inner = 2:(n-1)
    fixed_node = n
    c_fixed_val = 1.0
    
    F_eff = F[inner] - K[inner, fixed_node] .* c_fixed_val
    c_inner = K[inner, inner] \ Vector(F_eff)
    
    y_coeffs = zeros(n)
    y_coeffs[inner] = c_inner
    y_coeffs[n] = c_fixed_val
    
    return y_coeffs
end
