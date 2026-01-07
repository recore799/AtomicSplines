struct BSplineBasis
    order::Int
    num_splines::Int
    knots::Vector{Float64}
end

function bspline(i, k, t, knots)
    if k == 1
        return (knots[i] <= t < knots[i+1]) ? 1.0 : 0.0
    end
    d1 = knots[i+k-1] - knots[i]
    d2 = knots[i+k] - knots[i+1]
    v1 = (d1 == 0.0) ? 0.0 : (t - knots[i]) / d1 * bspline(i, k-1, t, knots)
    v2 = (d2 == 0.0) ? 0.0 : (knots[i+k] - t) / d2 * bspline(i+1, k-1, t, knots)
    return v1 + v2
end

function d_bspline(i, k, t, knots)
    if k == 1; return 0.0; end
    d1 = knots[i+k-1] - knots[i]
    d2 = knots[i+k] - knots[i+1]
    v1 = (d1 == 0.0) ? 0.0 : bspline(i, k-1, t, knots) / d1
    v2 = (d2 == 0.0) ? 0.0 : bspline(i+1, k-1, t, knots) / d2
    return (k-1) * (v1 - v2)
end

function eval_expansion(coeffs, basis, x)
    val = 0.0
    for i in 1:basis.num_splines
        if basis.knots[i] <= x < basis.knots[i+basis.order]
            val += coeffs[i] * bspline(i, basis.order, x, basis.knots)
        end
    end
    return val
end
