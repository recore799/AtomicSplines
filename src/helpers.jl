function evaluate_orbital(basis::BSplineBasis{K}, coeffs::AbstractVector{Float64}, r_points::AbstractVector{Float64}) where {K}
    n_points = length(r_points)
    psi_vals = zeros(Float64, n_points)
    
    vals = zeros(Float64, K)
    derivs = zeros(Float64, K)

    min_i = K
    max_i = length(basis.knots) - K + 1

    for (pt_idx, r) in enumerate(r_points)
        i = searchsortedlast(basis.knots, r)
        
        # Clamp to valid spline interval
        if i < min_i
            i = min_i
        elseif i > max_i
            i = max_i
        end

        eval_bspline_kernel!(vals, derivs, Val(true), Val(false), i, r, basis.knots, Val(K))
        
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
        push!(get!(shells, orb.l, Orbital[]), orb)
    end
    
    for (l, orbs) in shells
        sort!(orbs, by = x -> x.n)
    end
    
    return shells
end

