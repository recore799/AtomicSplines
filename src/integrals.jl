function init_scf_workspace(basis::BSplineBasis{K}, Z::Float64) where {K}
    n = basis.num_splines

    # Assemble core matrices
    S, T, V, V2, tensors = assemble_geometry(basis, Z)
    
    J = zeros(Float64, n, n)

    # Pre-allocate K matrices for s (l=0) and p (l=1) blocks.
    K_mats = Dict{Int, Matrix{Float64}}(
        0 => zeros(Float64, n, n),
        1 => zeros(Float64, n, n),
        2 => zeros(Float64, n, n),
        3 => zeros(Float64, n, n)
    )

    F_s = zeros(Float64, n,n)
    F_p = zeros(Float64, n,n)
    F_d = zeros(Float64, n,n)
    F_f = zeros(Float64, n,n)

    factors = Dict{Int, Any}()
    
    # Pre-allocate the N-length Poisson solver buffers
    b_buf = zeros(Float64, n)
    y_buf = zeros(Float64, n)
    y_orb_buffer = zeros(Float64, n)
    y_total_buffer = zeros(Float64, n)

    return SolverWorkspace{K}(
        basis, S, T, V, V2, J, 
        K_mats, F_s, F_p, F_d, F_f, tensors, factors, 
        zeros(Float64, K), zeros(Float64, K),
        b_buf, y_buf, y_orb_buffer, y_total_buffer
    )
end

function assemble_geometry(basis::BSplineBasis{K}, Z::Float64) where {K}
    n = basis.num_splines

    S  = zeros(Float64, n, n)
    T  = zeros(Float64, n, n)
    V  = zeros(Float64, n, n)
    V2 = zeros(Float64, n, n)   

    # Pre-allocate tensor storage
    tensors = Vector{Array{Float64, 3}}(undef, length(basis.knots) - 1)
    
    # Scratch buffers
    vals = zeros(Float64, K)
    derivs = zeros(Float64, K)

    # The Great Loop (Iterate Once)
    for i in 1:(length(basis.knots) - 1)
        t_a = basis.knots[i]; t_b = basis.knots[i+1]
        if t_a == t_b; continue; end
        
        mid = (t_a + t_b)/2; scale = (t_b - t_a)/2
        first_global = i - K + 1

        # Allocate local tensor for this element
        W_local = zeros(Float64, K, K, K)

        for q in 1:length(basis.gl_nodes)
            r = scale * basis.gl_nodes[q] + mid
            w = scale * basis.gl_weights[q]
            inv_r = (r > 1e-12) ? 1.0/r : 0.0
            inv_r2 = inv_r * inv_r

            # --- SINGLE KERNEL CALL ---
            eval_bspline_kernel!(vals, derivs, Val(true), Val(true), 
                                 i, r, basis.knots, Val(K))

            # --- ACCUMULATE EVERYTHING ---
            for a in 1:K
                g_a = first_global + a - 1
                
                # Pre-fetch values for 'a'
                Na = vals[a]; dNa = derivs[a]
                
                # Pre-calc parts of Tensor W to save mults
                # W_kab = (B_k/r) * B_a * B_b
                # We can cache (B_a * w * inv_r) for the tensor loops
                wa_tensor_factor = Na * w * inv_r 

                # Inner Loop (Triangular for Symmetry)
                for b in a:K
                    g_b = first_global + b - 1
                    
                    Nb = vals[b]; dNb = derivs[b]
                    
                    # Build Matrices S, T, V
                    # Check bounds only if strictly necessary
                    if g_a >= 1 && g_a <= n && g_b >= 1 && g_b <= n
                        term_S = w * Na * Nb
                        term_T = w * 0.5 * dNa * dNb
                        term_V = -Z * term_S * inv_r # Reuse term_S * inv_r
                        term_V2 = term_S * inv_r2
                        
                        S[g_a, g_b] += term_S
                        T[g_a, g_b] += term_T
                        V[g_a, g_b] += term_V
                        V2[g_a, g_b] += term_V2
                        
                    end
                    
                    # Build Tensor W (Using the same vals!)
                    # We need to loop 'k' here.
                    # W_kab = B_k * (B_a * B_b / r)
                    # We already have (B_a * B_b / r) implied.
                    term_tensor_base = wa_tensor_factor * Nb
                    
                    for k in 1:K
                        # W[k, a, b] += B_k * term_base
                        val = vals[k] * term_tensor_base
                        W_local[k, a, b] += val
                        
                        # Symmetry for Tensor W (a,b)
                        if a != b
                            W_local[k, b, a] += val
                        end
                    end
                end
            end
        end
        tensors[i] = W_local
    end
    
    return Symmetric(S), Symmetric(T), Symmetric(V), Symmetric(V2), tensors
end

function assemble_J_matrix!(ws::SolverWorkspace{K}, y_coeffs) where {K}
    n = ws.basis.num_splines
    J_mat = ws.J

    fill!(J_mat, 0.0)

    for i in 1:(length(ws.basis.knots) - 1)
        if !isassigned(ws.interaction_tensors, i); continue; end
        W_local = ws.interaction_tensors[i]
        
        first_global = i - K + 1
        
        y_local = MVector{K, Float64}(undef)
        for idx in 1:K
            g_idx = first_global + idx - 1
            y_local[idx] = (g_idx >= 1 && g_idx <= n) ? y_coeffs[g_idx] : 0.0
        end
        
        # Contract Tensor: J_ab = sum_k ( y_k * W_kab )
        for k_idx in 1:K
            weight = y_local[k_idx]
            if abs(weight) < 1e-12; continue; end
            
            for a in 1:K
                g_a = first_global + a - 1
                if g_a < 1 || g_a > n; continue; end
                
                for b_idx in 1:K
                    g_b = first_global + b_idx - 1
                    if g_b < 1 || g_b > n; continue; end

                    val_J = weight * W_local[k_idx, a, b_idx]
                    J_mat[g_a, g_b] += val_J
                end
            end
        end
    end
    
    return J_mat
end

function build_total_J_matrix!(ws::SolverWorkspace, orbitals::Vector{Orbital})
    fill!(ws.y_total_buffer, 0.0)
    for orb in orbitals
        if orb.occ > 0.0
            solve_poisson_J!(ws, ws.y_orb_buffer, orb)
            ws.y_total_buffer .+= orb.occ .* ws.y_orb_buffer
        end
    end
    
    return assemble_J_matrix!(ws, ws.y_total_buffer)
end

function assemble_K_matrix!(ws::SolverWorkspace{K}, K_mat::Matrix{Float64}, target_l::Int, orbitals::Vector{Orbital}) where {K}
    fill!(K_mat, 0.0)
    n = ws.basis.num_splines
    
    b_vec = zeros(Float64, n)
    y_k   = zeros(Float64, n)
    
    for psi in orbitals
        if psi.occ == 0.0; continue; end
        
        # --- Determine Allowed Multipoles & Multipliers ---
        multipoles = Tuple{Int, Float64}[]
        if target_l == 0 && psi.l == 0
            push!(multipoles, (0, 1.0))
        elseif target_l == 0 && psi.l == 1
            push!(multipoles, (1, 1.0 / 3.0))
        elseif target_l == 1 && psi.l == 0
            push!(multipoles, (1, 1.0 / 3.0))
        elseif target_l == 1 && psi.l == 1
            push!(multipoles, (0, 1.0 / 3.0))
            push!(multipoles, (2, 2.0 / 15.0))
        # --- s-d and d-s ---
        elseif (target_l == 0 && psi.l == 2) || (target_l == 2 && psi.l == 0)
            push!(multipoles, (2, 1.0 / 5.0))
            
        # --- p-d and d-p ---
        elseif (target_l == 1 && psi.l == 2) || (target_l == 2 && psi.l == 1)
            push!(multipoles, (1, 2.0 / 15.0))
            push!(multipoles, (3, 3.0 / 35.0))
            
        # --- d-d ---
        elseif target_l == 2 && psi.l == 2
            push!(multipoles, (0, 1.0 / 5.0))
            push!(multipoles, (2, 2.0 / 35.0))
            push!(multipoles, (4, 2.0 / 35.0))
        
        # --- s-f and f-s ---
        elseif (target_l == 0 && psi.l == 3) || (target_l == 3 && psi.l == 0)
            push!(multipoles, (3, 1.0 / 7.0))
            
        # --- p-f and f-p ---
        elseif (target_l == 1 && psi.l == 3) || (target_l == 3 && psi.l == 1)
            push!(multipoles, (2, 3.0 / 35.0))
            push!(multipoles, (4, 4.0 / 63.0))
            
        # --- d-f and f-d ---
        elseif (target_l == 2 && psi.l == 3) || (target_l == 3 && psi.l == 2)
            push!(multipoles, (1, 3.0 / 35.0))
            push!(multipoles, (3, 4.0 / 105.0))
            push!(multipoles, (5, 10.0 / 231.0))
            
        # --- f-f ---
        elseif target_l == 3 && psi.l == 3
            push!(multipoles, (0, 1.0 / 7.0))
            push!(multipoles, (2, 4.0 / 105.0))
            push!(multipoles, (4, 2.0 / 77.0))
            push!(multipoles, (6, 100.0 / 3003.0))
        end
        
        spin_weight = psi.occ / 2.0
        
        for nu in 1:n
            fill!(b_vec, 0.0)
            
            # --- Build Source Vector ---
            for i in 1:(length(ws.basis.knots)-1)
                if !isassigned(ws.interaction_tensors, i); continue; end
                W = ws.interaction_tensors[i]
                first = i - K + 1
                
                loc_nu = nu - first + 1
                if loc_nu < 1 || loc_nu > K; continue; end
                
                # Zero-allocation stack vector
                c_loc = MVector{K, Float64}(undef)
                for idx in 1:K
                    g = first + idx - 1
                    c_loc[idx] = (g >= 1 && g <= n) ? psi.coeffs[g] : 0.0
                end
                
                for k_idx in 1:K
                    g_k = first + k_idx - 1
                    if g_k < 1 || g_k > n; continue; end
                    val = 0.0
                    for a in 1:K
                        val += c_loc[a] * W[k_idx, a, loc_nu]
                    end
                    b_vec[g_k] += val
                end
            end
            
            # --- Solve & Accumulate ---
            for (k_mult, angular_mult) in multipoles
                # Mutates y_k directly!
                solve_generalized_poisson!(ws, y_k, b_vec, k_mult)
                total_multiplier = angular_mult * spin_weight
                
                for i in 1:(length(ws.basis.knots)-1)
                    if !isassigned(ws.interaction_tensors, i); continue; end
                    W = ws.interaction_tensors[i]
                    first = i - K + 1

                    # Zero-allocation stack vectors!
                    c_psi = MVector{K, Float64}(undef)
                    c_y   = MVector{K, Float64}(undef)
                    
                    for idx in 1:K
                        g = first + idx - 1
                        if g >= 1 && g <= n
                            c_psi[idx] = psi.coeffs[g]
                            c_y[idx]   = y_k[g]
                        else
                            c_psi[idx] = 0.0
                            c_y[idx]   = 0.0
                        end
                    end

                    for mu in 1:K
                        g_mu = first + mu - 1
                        if g_mu < 1 || g_mu > n; continue; end
                        
                        val = 0.0
                        for a in 1:K
                            pot = 0.0
                            for lam in 1:K
                                pot += c_y[lam] * W[lam, mu, a]
                            end
                            val += c_psi[a] * pot
                        end
                        K_mat[g_mu, nu] += val * total_multiplier 
                    end
                end
            end
        end 
    end
    
    # for i in 1:n
    #     for j in (i+1):n
    #         K_mat[j, i] = K_mat[i, j]
    #     end
    # end
    
    return Symmetric(K_mat)
end
