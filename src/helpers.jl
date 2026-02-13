function evaluate_orbital(basis, coeffs, r_points::AbstractVector{Float64})
    n_points = length(r_points)
    psi_vals = zeros(Float64, n_points)
    k = basis.order
    
    # Pre-allocate scratch
    scratch_vals = zeros(Float64, k)
    scratch_derivs = zeros(Float64, k)

    current_element = 1
    
    # CRITICAL FIX: The last element we can evaluate is limited by the spline order.
    # If we go past this, the kernel will try to read knots that don't exist.
    last_safe_element = length(basis.knots) - k

    for (pt_idx, r) in enumerate(r_points)
        
        # Reset if the grid is not sorted (safety check)
        if r < basis.knots[current_element]
            current_element = 1
        end

        # 1. Find interval. 
        # FIX: STRICT upper bound check (last_safe_element)
        while current_element < last_safe_element && r >= basis.knots[current_element+1]
            current_element += 1
        end
        
        # 2. Identify 'i'
        i = current_element
        
        # 3. Kernel Call
        # If r > R_max, this might extrapolate slightly or give 0, but it won't crash.
        eval_bspline_kernel!(scratch_vals, scratch_derivs, Val(true), Val(false), 
                             i, k, r, basis.knots)
        
        # 4. Dot product with coefficients
        first_global = i - k + 1
        val = 0.0
        for local_idx in 1:k
            global_idx = first_global + local_idx - 1
            if global_idx >= 1 && global_idx <= basis.num_splines
                val += coeffs[global_idx] * scratch_vals[local_idx]
            end
        end
        psi_vals[pt_idx] = val
    end
    
    return psi_vals
end
