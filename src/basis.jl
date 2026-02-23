# K is the spline order. It is now a compile-time constant.
struct BSplineBasis{K}
    num_splines::Int
    knots::Vector{Float64}
    
    # Pre-computed Gauss-Legendre nodes/weights for this order
    # We store them here to avoid 'gausslegendre(k+4)' calls inside loops
    gl_nodes::Vector{Float64}
    gl_weights::Vector{Float64}
end

# Holds all the heavy allocations so the solver loop is allocation-free.
struct SolverWorkspace{K}
    basis::BSplineBasis{K}
    
    # S, T, V have bandwidth K-1 each side.
    S::BandedMatrix{Float64}
    T::BandedMatrix{Float64}
    V::BandedMatrix{Float64}
    V2::BandedMatrix{Float64}
    
    # Coulomb matrix
    J::BandedMatrix{Float64}
    
    # Maps angular momentum 'l' to its specific K matrix
    K_mats::Dict{Int, Matrix{Float64}}

    # Stores the local KxKxK tensor for each element
    interaction_tensors::Vector{Array{Float64, 3}}

    # Store multiple Cholesky factorizations
    poisson_factors::Dict{Int, Any}

    # Pre-allocated Scratch Vectors (Size K)
    scratch_vals::Vector{Float64}
    scratch_derivs::Vector{Float64}
    
    # Size N buffers for Poisson solvers ---
    b_buffer::Vector{Float64}
    y_buffer::Vector{Float64}
end

# struct SolverWorkspace{K}
#     basis::BSplineBasis{K}
    
#     # S, T, V have bandwith K-1 each side.
#     S::BandedMatrix{Float64}
#     T::BandedMatrix{Float64}
#     V::BandedMatrix{Float64}
#     V2::BandedMatrix{Float64}
#     # Coulomb matrix
#     J::BandedMatrix{Float64}
#     K_mat::Matrix{Float64}

   
#     # Stores the local KxKxK tensor for each element
#     # W[element_index][k, a, b]
#     interaction_tensors::Vector{Array{Float64, 3}}

#     # Store multiple Cholesky factorizations
#     # Key = k (multipole order), Value = Factorization
#     poisson_factors::Dict{Int, Any}

#     # Pre-allocated Scratch Vectors (Size K)
#     scratch_vals::Vector{Float64}
#     scratch_derivs::Vector{Float64}
    
# end

function generate_basis(R_max, N_elems, ::Val{K}; γ=2.0) where {K}
    # Knot Generation
    x = range(0, 1, length=N_elems + 1)
    r_knots = R_max .* (exp.(γ .* x) .- 1) ./ (exp(γ) - 1)
    r_knots[1] = 0.0 # Clamp small errors to 0
    
    # Open knot vector (repeat start/end K-1 times)
    knots = vcat(fill(r_knots[1], K-1), r_knots, fill(r_knots[end], K-1))
    
    n_splines = length(knots) - K
    
    # 2. Pre-compute Quadrature (Order + 6 rule)
    # This runs once at setup, so allocation is fine here.
    gl_p, gl_w = gausslegendre(K + 6)
    
    return BSplineBasis{K}(n_splines, knots, gl_p, gl_w)
end

"""
Corrected Parametric Kernel.
Calculates derivatives using Order K-1 values BEFORE overwriting them with Order K.
"""
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

    # Initialize Order 1
    @inbounds out_vals[1] = 1.0

    # Cox-de Boor Loop: Run ONLY up to Order K-1
    # We stop at j = K-2 so out_vals holds B_{K-1}
    for j in 1:(K-2)
        saved = 0.0
        for r in 1:j
            left  = knots[i - j + r]
            right = knots[i + r]
            if right == left
                term = 0.0
            else
                term = out_vals[r] / (right - left)
            end
            out_vals[r] = saved + (right - t) * term
            saved = (t - left) * term
        end
        out_vals[j+1] = saved
    end

    if CompDeriv
        # We are computing B'_{K}. This depends on B_{K-1} (currently in out_vals)
        for r in 1:K
            # Handle boundary logic for the indices
            val_im1 = (r == 1) ? 0.0 : out_vals[r-1]
            val_i   = (r == K) ? 0.0 : out_vals[r]

            denom1 = knots[i+r-1] - knots[i-K+r]
            denom2 = knots[i+r]   - knots[i-K+r+1]

            term1 = (denom1 == 0.0) ? 0.0 : val_im1 / denom1
            term2 = (denom2 == 0.0) ? 0.0 : val_i   / denom2

            out_derivs[r] = (K - 1) * (term1 - term2)
        end
    end

    # Finish Computing Final Values (Order K)
    # Now safe to overwrite out_vals with the final step
    if CompVal
        j = K - 1
        saved = 0.0
        for r in 1:j
            left  = knots[i - j + r]
            right = knots[i + r]
            if right == left
                term = 0.0
            else
                term = out_vals[r] / (right - left)
            end
            out_vals[r] = saved + (right - t) * term
            saved = (t - left) * term
        end
        out_vals[j+1] = saved
    end
end
