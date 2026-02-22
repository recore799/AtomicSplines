function init_scf_workspace(basis::BSplineBasis{K}, Z::Float64) where {K}
    n = basis.num_splines

    # println("   > Assembling Geometry (Matrices & Tensors)...")
    S, T, V, V2, tensors = assemble_geometry(basis, Z)
    
    bw = (K-1, K-1)
    J = BandedMatrix(Zeros(n, n), bw)

    # Initialize Factorization Dictionary
    factors = Dict{Int, Any}()
    
    # Pre-compute k=0 (Standard Poisson)
    #    Operator: 2*T + 0*V2
    # println("   > Factorizing Poisson (k=0)...")
    K0 = 2.0 .* T
    active = 2:(n-1)
    factors[0] = cholesky(K0[active, active])

    K_mat = zeros(n,n)
   
    return SolverWorkspace(basis, S, T, V, V2, J, K_mat, tensors, 
                           factors, zeros(K), zeros(K))
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

function assemble_K_matrix(ws::SolverWorkspace{K}, target_l::Int, orbitals::Vector{Orbital}) where {K}
    K_mat = ws.K_mat
    fill!(K_mat, 0.0)
    n = ws.basis.num_splines
    b_vec = zeros(Float64, n)
    
    for psi in orbitals
        if psi.occ == 0.0; continue; end
        
        # --- Determine Allowed Multipoles & Multipliers ---
        # Tuples of (k_mult, angular_multiplier)
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
        
        # Exchange weight accounts for same-spin interactions (occ / 2.0)
        spin_weight = psi.occ / 2.0
        
        for nu in 1:n
            fill!(b_vec, 0.0)
            
            # --- Compute Boundary Condition (Q) ---
            # Used only for k=0 evaluations
            boundary_charge = 0.0
            s_start = max(1, nu - K + 1)
            s_end   = min(n, nu + K - 1)
            for a in s_start:s_end
                boundary_charge += psi.coeffs[a] * ws.S[a, nu]
            end
            
            # --- Build Source Vector ---
            for i in 1:(length(ws.basis.knots)-1)
                if !isassigned(ws.interaction_tensors, i); continue; end
                W = ws.interaction_tensors[i]
                first = i - K + 1
                
                loc_nu = nu - first + 1
                if loc_nu < 1 || loc_nu > K; continue; end
                
                c_loc = zeros(Float64, K)
                for idx in 1:K
                    g = first + idx - 1
                    if g>=1 && g<=n; c_loc[idx] = psi.coeffs[g]; end
                end
                
                for k_idx in 1:K
                    g_k = first + k_idx - 1; if g_k<1 || g_k>n; continue; end
                    val = 0.0
                    for a in 1:K
                        val += c_loc[a] * W[k_idx, a, loc_nu]
                    end
                    b_vec[g_k] += val
                end
            end
            
            # --- Solve & Accumulate for Each Multipole ---
            for (k_mult, angular_mult) in multipoles
                
                y_k = solve_generalized_poisson(ws, b_vec, k_mult)
                total_multiplier = angular_mult * spin_weight
                
                # --- Accumulate into K_mat ---
                for i in 1:(length(ws.basis.knots)-1)
                     if !isassigned(ws.interaction_tensors, i); continue; end
                     W = ws.interaction_tensors[i]
                     first = i - K + 1

                     c_psi = zeros(Float64, K); c_y = zeros(Float64, K)
                     for idx in 1:K
                         g = first + idx - 1
                         if g>=1 && g<=n
                             c_psi[idx] = psi.coeffs[g]
                             c_y[idx] = y_k[g]
                         end
                     end

                     for mu in 1:K
                         g_mu = first + mu - 1
                         if g_mu < 1 || g_mu > n; continue; end
                         
                         val = 0.0
                         for a in 1:K
                             pot = 0.0
                             for lam in 1:K; pot += c_y[lam] * W[lam, mu, a]; end
                             val += c_psi[a] * pot
                         end
                         
                         # Apply angular and spin weights to the exchange integral
                         K_mat[g_mu, nu] += val * total_multiplier 
                     end
                end
            end
        end 
    end
    return Symmetric(copy(K_mat))
end
