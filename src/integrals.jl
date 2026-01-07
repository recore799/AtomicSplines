function assemble_core(basis, Z)
    n = basis.num_splines; k = basis.order; knots = basis.knots
    T = spzeros(n, n); V_nuc = spzeros(n, n); S = spzeros(n, n)
    gl_p, gl_w = gausslegendre(k + 2)

    for i in 1:n
        for j in i:min(i+k-1, n)
            t_val = 0.0; v_val = 0.0; s_val = 0.0
            
            s_knot = max(i, j); e_knot = min(i+k, j+k)
            for interval in s_knot:(e_knot-1)
                ta = knots[interval]; tb = knots[interval+1]; if ta==tb continue end
                mid = (ta+tb)/2; scale = (tb-ta)/2
                for q in 1:length(gl_p)
                    r = scale*gl_p[q] + mid; w = scale*gl_w[q]
                    
                    bi = bspline(i, k, r, knots); bj = bspline(j, k, r, knots)
                    dbi = d_bspline(i, k, r, knots); dbj = d_bspline(j, k, r, knots)
                    
                    s_val += w * bi * bj
                    t_val += w * 0.5 * dbi * dbj
                    v_val += w * (-Z/r) * bi * bj
                end
            end
            
            S[i,j] = S[j,i] = s_val
            T[i,j] = T[j,i] = t_val
            V_nuc[i,j] = V_nuc[j,i] = v_val
        end
    end
    return T, V_nuc, S
end

function assemble_J_matrix(basis, y_coeffs)
    n = basis.num_splines; k = basis.order
    J_mat = spzeros(n, n)
    gl_p, gl_w = gausslegendre(k + 2)

    for i in 1:n
        for j in i:min(i+k-1, n)
            val = 0.0
            s_knot = max(i, j); e_knot = min(i+k, j+k)
            for interval in s_knot:(e_knot-1)
                ta = basis.knots[interval]; tb = basis.knots[interval+1]; if ta==tb continue end
                mid = (ta+tb)/2; scale = (tb-ta)/2
                for q in 1:length(gl_p)
                    r = scale*gl_p[q] + mid; w = scale*gl_w[q]
                    y_val = eval_expansion(y_coeffs, basis, r)
                    V_H = y_val / r
                    val += w * V_H * bspline(i, k, r, basis.knots) * bspline(j, k, r, basis.knots)
                end
            end
            J_mat[i,j] = J_mat[j,i] = val
        end
    end
    return J_mat
end
