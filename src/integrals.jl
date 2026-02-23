function init_scf_workspace(basis::BSplineBasis{K}, Z::Float64) where {K}
    n = basis.num_splines

    # Assemble core matrices
    S, T, V, V2, tensors = assemble_geometry(basis, Z)
    
    bw = (K-1, K-1)
    J = BandedMatrix(Zeros(n, n), bw)

    # Pre-allocate K matrices for s (l=0) and p (l=1) blocks.
    K_mats = Dict{Int, Matrix{Float64}}(
        0 => zeros(Float64, n, n),
        1 => zeros(Float64, n, n)
    )

    factors = Dict{Int, Any}()
    
    # Pre-allocate the N-length Poisson solver buffers
    b_buf = zeros(Float64, n)
    y_buf = zeros(Float64, n)

    return SolverWorkspace{K}(
        basis, S, T, V, V2, J, 
        K_mats, tensors, factors, 
        zeros(Float64, K), zeros(Float64, K),
        b_buf, y_buf
    )
end

function assemble_geometry(basis::BSplineBasis{K}, Z::Float64) where {K}
    n = basis.num_splines
    bw = (K-1, K-1)
    
    #    We use 'Zeros' which BandedMatrices understands efficiently
    S = BandedMatrix(Zeros(n, n), bw)
    T = BandedMatrix(Zeros(n, n), bw)
    V = BandedMatrix(Zeros(n, n), bw)
    V2 = BandedMatrix(Zeros(n, n), bw)
    
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
                        
                        if g_a != g_b
                            S[g_b, g_a] += term_S
                            T[g_b, g_a] += term_T
                            V[g_b, g_a] += term_V
                            V2[g_b, g_a] += term_V2
                        end
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
    
    return S, T, V, V2, tensors
end

function assemble_J_matrix(ws::SolverWorkspace{K}, y_coeffs) where {K}
    n = ws.basis.num_splines
    J_mat = ws.J

    fill!(J_mat, 0.0)

    for i in 1:(length(ws.basis.knots) - 1)
        if !isassigned(ws.interaction_tensors, i); continue; end
        W_local = ws.interaction_tensors[i]
        
        first_global = i - K + 1
        
        # Extract local y coefficients
        y_local = zeros(Float64, K)
        for idx in 1:K
            g_idx = first_global + idx - 1
            if g_idx >= 1 && g_idx <= n
                y_local[idx] = y_coeffs[g_idx]
            end
        end
        
        # Contract Tensor: J_ab = sum_k ( y_k * W_kab )
        for k in 1:K
            weight = y_local[k]
            if abs(weight) < 1e-12; continue; end
            
            for a in 1:K
                g_a = first_global + a - 1
                if g_a < 1 || g_a > n; continue; end
                
                for b in 1:K
                    g_b = first_global + b - 1
                    if g_b < 1 || g_b > n; continue; end

                    val_J = weight * W_local[k, a, b]
                    J_mat[g_a, g_b] += val_J
                end
            end
        end
    end
    
    return J_mat
end

function build_total_J_matrix(ws::SolverWorkspace, orbitals::Vector{Orbital})
    y_total = zeros(Float64, ws.basis.num_splines)
    for orb in orbitals
        if orb.occ > 0.0
            y_orb = solve_poisson_J(ws, orb)
            y_total .+= orb.occ .* y_orb
        end
    end
    
    return assemble_J_matrix(ws, y_total)
end

using StaticArrays

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
    
    # Explicitly symmetrize in-place to avoid `copy()` and type changes
    for i in 1:n
        for j in (i+1):n
            K_mat[j, i] = K_mat[i, j]
        end
    end
    
    return K_mat
end
