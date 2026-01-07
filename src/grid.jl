function generate_basis(r_max, num_elements, order; γ=2.0)
    x = range(0, 1, length=num_elements + 1)
    r_knots = r_max .* (exp.(γ .* x) .- 1) ./ (exp(γ) - 1)
    knots = vcat(fill(r_knots[1], order-1), r_knots, fill(r_knots[end], order-1))
    return BSplineBasis(order, length(knots)-order, knots)
end
