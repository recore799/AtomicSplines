function solve_generalized_poisson!(ws::SolverWorkspace{K}, y_full::Vector{Float64}, source::Vector{Float64}, k_mult::Int) where {K}
    start_idx = k_mult + 2
    n = ws.basis.num_splines
    active_k = start_idx:n
    length_active = length(active_k)
    
    # Lazy Factorization
    if !haskey(ws.poisson_factors, k_mult)
        A_k = 2 * ws.T .+ (k_mult * (k_mult + 1)) .* ws.V2
        A_k_dense = Matrix(A_k[active_k, active_k])
        A_k_dense[end, end] += k_mult / ws.basis.knots[end] # Robin BC (or is it Noumann?)
        ws.poisson_factors[k_mult] = cholesky(Symmetric(A_k_dense))
    end
    
    # View into workspace buffers to prevent slicing allocations
    b_view = view(ws.b_buffer, 1:length_active)
    y_view = view(ws.y_buffer, 1:length_active)
    
    # Scale source by (2k + 1) directly into the buffer
    for (i, g) in enumerate(active_k)
        b_view[i] = source[g] * (2.0 * k_mult + 1.0)
    end
    
    # In-place Cholesky solve: writes result directly into y_view
    ldiv!(y_view, ws.poisson_factors[k_mult], b_view)
    
    # Update y_full
    fill!(y_full, 0.0)
    for (i, g) in enumerate(active_k)
        y_full[g] = y_view[i]
    end
    
    return y_full
end


# function solve_generalized_poisson!(ws::SolverWorkspace{K}, y_full::Vector{Float64}, source::Vector{Float64}, k_mult::Int) where {K}
#     start_idx = k_mult + 2
#     n = ws.basis.num_splines
#     active_k = start_idx:n
    
#     # Lazy factorization
#     if !haskey(ws.poisson_factors, k_mult)
#         A_k = 2 * ws.T .+ (k_mult * (k_mult + 1)) .* ws.V2
#         A_k_dense = Matrix(A_k[active_k, active_k])
        
#         # Apply Robin (or is it Neuman?) Boundary Condition at the wall
#         A_k_dense[end, end] += k_mult / R_max
        
#         ws.poisson_factors[k_mult] = cholesky(Symmetric(A_k_dense))
#     end
    
#     # non-allocation
#     length_active = length(active_k)
#     b_view = view(ws.b_buffer, 1:length_active)
    
#     for (i, g) in enumerate(active_k)
#         b_view[i] = source[g] * (2.0 * k_mult + 1)
#     end
    
#     # In-place linear solve
#     y_view = view(ws.y_buffer, 1:length_active)
    
#     # ldiv! writes the result of (A \ b) directly into y_view
#     ldiv!(y_view, ws.poisson_factors[k_mult], b_view)
    
#     # Write back to the mutated output array
#     fill!(y_full, 0.0)
#     y_full[active_k] .= y_view
    
#     return y_full
# end
# function solve_generalized_poisson(ws::SolverWorkspace{K}, source::Vector{Float64}, k_mult::Int) where {K}
#     start_idx = k_mult + 2
#     n = ws.basis.num_splines
#     R_max = ws.basis.knots[end]
    
#     # Extend active set to the boundary 'n'
#     active_k = start_idx:n
    
#     if !haskey(ws.poisson_factors, k_mult)
#         A_k = 2 * ws.T .+ (k_mult * (k_mult + 1)) .* ws.V2
#         A_k_dense = Matrix(A_k[active_k, active_k])
        
#         # Apply Robin (or is it Neuman?) Boundary Condition at the wall
#         A_k_dense[end, end] += k_mult / R_max
        
#         ws.poisson_factors[k_mult] = cholesky(Symmetric(A_k_dense))
#     end
    
#     b_inner = source[active_k] .* (2.0 * k_mult + 1)
#     y_inner = ws.poisson_factors[k_mult] \ b_inner
    
#     y_full = zeros(Float64, n)
#     y_full[active_k] = y_inner
    
#     return y_full
# end

function solve_poisson_J(ws::SolverWorkspace{K}, orb::Orbital) where {K}
    # Coulomb potentials for closed shells are spherically symmetric (k=0)
    basis = ws.basis; n = basis.num_splines
    b = zeros(Float64, n); c_local = ws.scratch_vals
    
    for i in 1:(length(basis.knots)-1)
        if !isassigned(ws.interaction_tensors, i); continue; end
        W_local = ws.interaction_tensors[i]
        first_global = i - K + 1
        
        fill!(c_local, 0.0)
        for idx in 1:K
            g = first_global + idx - 1
            if g >= 1 && g <= n; c_local[idx] = orb.coeffs[g]; end
        end
        
        for k_idx in 1:K
            g_k = first_global + k_idx - 1
            if g_k < 1 || g_k > n; continue; end
            val = 0.0
            for a in 1:K
                ca = c_local[a]
                if abs(ca) < 1e-15; continue; end
                val += ca * ca * W_local[k_idx, a, a]
                for b in (a+1):K
                    val += 2.0 * ca * c_local[b] * W_local[k_idx, a, b]
                end
            end
            b[g_k] += val
        end
    end
    
    active = 2:n
    
    if !haskey(ws.poisson_factors, 0)
        A_0_dense = 2 .* Matrix(ws.T)[active, active]
        ws.poisson_factors[0] = cholesky(Symmetric(A_0_dense))
    end
    
    b_inner = b[active] 
    y_inner = ws.poisson_factors[0] \ b_inner
    
    y_full = zeros(Float64, n)
    y_full[active] = y_inner
    
    return y_full
end
