function evaluate_orbital(basis::BSplineBasis{K}, coeffs::AbstractVector{Float64}, r_points::AbstractVector{Float64}) where {K}
    n_points = length(r_points)
    psi_vals = zeros(Float64, n_points)
    
    # Scratch buffers
    vals = zeros(Float64, K)
    derivs = zeros(Float64, K)

    # Calculate strict bounds for the interval index 'i'
    # The kernel at 'i' accesses knots up to 'i + K - 1'
    # So we need: i + K - 1 <= length(knots)  =>  i <= length - K + 1
    min_i = K
    max_i = length(basis.knots) - K + 1

    for (pt_idx, r) in enumerate(r_points)
        
        # 1. Find Knot Interval
        # Returns index 'i' such that knots[i] <= r < knots[i+1]
        i = searchsortedlast(basis.knots, r)
        
        # 2. CLAMPING (The Fix)
        if i < min_i
            i = min_i
        elseif i > max_i
            i = max_i
        end

        # 3. Kernel Call
        eval_bspline_kernel!(vals, derivs, Val(true), Val(false), 
                             i, r, basis.knots, Val(K))
        
        # 4. Contract with Coefficients
        # The splines active in interval 'i' are indices [i-K+1 ... i]
        first_global = i - K + 1
        val = 0.0
        
        for local_idx in 1:K
            global_idx = first_global + local_idx - 1
            if global_idx >= 1 && global_idx <= basis.num_splines
                val += coeffs[global_idx] * vals[local_idx]
            end
        end
        psi_vals[pt_idx] = val
    end
    
    return psi_vals
end

function group_orbitals_by_l(orbitals::Vector{Orbital})
    shells = Dict{Int, Vector{Orbital}}()
    for orb in orbitals
        # If the l-key doesn't exist, create an empty array and push the orbital
        push!(get!(shells, orb.l, Orbital[]), orb)
    end
    
    # Sort them by principal quantum number 'n' within each l-block
    for (l, orbs) in shells
        sort!(orbs, by = x -> x.n)
    end
    
    return shells
end
