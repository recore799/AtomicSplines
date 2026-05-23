struct BSplineBasis{K}
    num_splines::Int
    knots::Vector{Float64}
    gl_nodes::Vector{Float64}
    gl_weights::Vector{Float64}
end

struct SolverWorkspace{K}
    basis::BSplineBasis{K}
    
    S::Matrix{Float64}
    T::Matrix{Float64}
    V::Matrix{Float64}
    R::Matrix{Float64}
    R2::Matrix{Float64}
    R_inv2::Matrix{Float64}
    R_inv3::Matrix{Float64}
    
    J::Matrix{Float64}
    
    K_mats::Dict{Int, Matrix{Float64}}
    F_mats::Dict{Int, Matrix{Float64}}

    interaction_tensors::Vector{Array{Float64, 3}}
    poisson_factors::Dict{Int, Any}

    scratch_vals::Vector{Float64}
    scratch_derivs::Vector{Float64}
    
    b_buffer::Vector{Float64}
    y_buffer::Vector{Float64}
    y_orb_buffer::Vector{Float64}   
    y_total_buffer::Vector{Float64} 

    scratch_source::Vector{Float64}
    scratch_y::Vector{Float64}
    rk_cache::Dict{UInt64, Float64}
end

function generate_basis(R_max, N_elems, ::Val{K}; γ=2.0) where {K}
    x = range(0, 1, length=N_elems + 1)
    r_knots = R_max .* (exp.(γ .* x) .- 1) ./ (exp(γ) - 1)
    r_knots[1] = 0.0 
    
    knots = vcat(fill(r_knots[1], K-1), r_knots, fill(r_knots[end], K-1))
    n_splines = length(knots) - K
    
    gl_p, gl_w = gausslegendre(K + 6)
    
    return BSplineBasis{K}(n_splines, knots, gl_p, gl_w)
end

@inline function eval_bspline_kernel!(
    out_vals::AbstractVector{Float64}, 
    out_derivs::AbstractVector{Float64},
    ::Val{CompVal}, 
    ::Val{CompDeriv}, 
    i::Int, 
    t::Float64, 
    knots::Vector{Float64},
    ::Val{K}
) where {CompVal, CompDeriv, K}

    @inbounds out_vals[1] = 1.0

    for j in 1:(K-2)
        saved = 0.0
        for r in 1:j
            left  = knots[i - j + r]
            right = knots[i + r]
            term = (right == left) ? 0.0 : out_vals[r] / (right - left)
            out_vals[r] = saved + (right - t) * term
            saved = (t - left) * term
        end
        out_vals[j+1] = saved
    end

    if CompDeriv
        for r in 1:K
            val_im1 = (r == 1) ? 0.0 : out_vals[r-1]
            val_i   = (r == K) ? 0.0 : out_vals[r]

            denom1 = knots[i+r-1] - knots[i-K+r]
            denom2 = knots[i+r]   - knots[i-K+r+1]

            term1 = (denom1 == 0.0) ? 0.0 : val_im1 / denom1
            term2 = (denom2 == 0.0) ? 0.0 : val_i   / denom2

            out_derivs[r] = (K - 1) * (term1 - term2)
        end
    end

    if CompVal
        j = K - 1
        saved = 0.0
        for r in 1:j
            left  = knots[i - j + r]
            right = knots[i + r]
            term = (right == left) ? 0.0 : out_vals[r] / (right - left)
            out_vals[r] = saved + (right - t) * term
            saved = (t - left) * term
        end
        out_vals[j+1] = saved
    end
end
