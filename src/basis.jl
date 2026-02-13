struct BSplineBasis
    order::Int
    num_splines::Int
    knots::Vector{Float64}
    num_elements::Int
end

struct BSplineScratch
    vals::Vector{Float64}
    derivs::Vector{Float64}
end

"""
Unified B-spline evaluation kernel.
Computes Values (if CompVal) and Derivatives (if CompDeriv) in one pass.
Writes directly to pre-allocated scratch buffers to avoid GC.
"""
function eval_bspline_kernel!(
    out_vals::AbstractVector{Float64}, 
    out_derivs::AbstractVector{Float64},
    ::Val{CompVal}, 
    ::Val{CompDeriv}, 
    i::Int, k::Int, t::Float64, knots::Vector{Float64}
) where {CompVal, CompDeriv}

    # 1. Initialize Order 1
    out_vals[1] = 1.0

    # 2. Cox-de Boor recursion up to Order k-1
    for j in 1:(k-2)
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

    # 3. Compute Derivatives (Optional) - uses Order k-1 values
    if CompDeriv
        for r in 1:k
            val_im1 = (r == 1) ? 0.0 : out_vals[r-1]
            val_i   = (r == k) ? 0.0 : out_vals[r]

            denom1 = knots[i+r-1] - knots[i-k+r]
            denom2 = knots[i+r]   - knots[i-k+r+1]

            term1 = (denom1 == 0.0) ? 0.0 : val_im1 / denom1
            term2 = (denom2 == 0.0) ? 0.0 : val_i   / denom2

            out_derivs[r] = (k - 1) * (term1 - term2)
        end
    end

    # 4. Compute Final Values (Optional) - Order k-1 -> k
    if CompVal
        j = k - 1
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

function generate_basis(r_max, num_elements, order; γ=2.0)
    x = range(0, 1, length=num_elements + 1)
    r_knots = r_max .* (exp.(γ .* x) .- 1) ./ (exp(γ) - 1)
    # Clamp small errors to 0
    r_knots[1] = 0.0
    
    # Open knot vector (repeat start and end)
    knots = vcat(fill(r_knots[1], order-1), r_knots, fill(r_knots[end], order-1))
    return BSplineBasis(order, length(knots)-order, knots, length(knots)-1)
end

