function init_scf_workspace(basis::BSplineBasis{K}, Z::Float64; calc_R_matrices::Bool=false) where {K}
    n = basis.num_splines

    S, T, V, R_inv2, R_inv3, R, R2, tensors = assemble_geometry(basis, Z; calc_R_matrices=calc_R_matrices)
    
    J = zeros(Float64, n, n)

    K_mats = Dict{Int, Matrix{Float64}}(l => zeros(Float64, n, n) for l in 0:5)
    F_mats = Dict{Int, Matrix{Float64}}(l => zeros(Float64, n, n) for l in 0:5)

    factors = Dict{Int, Any}()
    
    b_buf = zeros(Float64, n)
    y_buf = zeros(Float64, n)
    y_orb_buffer = zeros(Float64, n)
    y_total_buffer = zeros(Float64, n)

    scratch_source = zeros(Float64, n)
    scratch_y = zeros(Float64, n)
    rk_cache = Dict{UInt64, Float64}()

    return SolverWorkspace{K}(
        basis, S, T, V, R_inv2, R_inv3, R, R2, J, 
        K_mats, F_mats, tensors, factors, 
        zeros(Float64, K), zeros(Float64, K),
        b_buf, y_buf, y_orb_buffer, y_total_buffer, scratch_source, scratch_y, rk_cache
    )
end


function cached_init_scf_workspace(R_max::Float64, N_elems::Int, ::Val{K}, Z::Float64; γ::Float64=2.0, calc_R_matrices::Bool=true) where {K}

    filename = @sprintf("geometry_Z%.1f_R%.1f_N%d_K%d_g%.2f.jld2", Z, R_max, N_elems, K, γ)
    
    basis = generate_basis(R_max, N_elems, Val(K); γ=γ)
    n = basis.num_splines
    
    if isfile(filename)
        println("Loading cached geometry from $filename ...")
        data = load(filename)
        S = data["S"]
        T = data["T"]
        V = data["V"]
        R = data["R"]
        R2 = data["R2"]
        R_inv2 = data["R_inv2"]
        R_inv3 = data["R_inv3"]
        tensors = data["tensors"]
    else
        println("No cache found. Computing geometry and saving to $filename ...")
        S, T, V, R, R2, R_inv2, R_inv3, tensors = assemble_geometry(basis, Z; calc_R_matrices=calc_R_matrices)
        
        jldsave(filename; S=S, T=T, V=V, R=R, R2=R2, R_inv2=R_inv2, R_inv3=R_inv3, tensors=tensors)
    end
    
    J = zeros(Float64, n, n)
    K_mats = Dict{Int, Matrix{Float64}}(l => zeros(Float64, n, n) for l in 0:5)
    F_mats = Dict{Int, Matrix{Float64}}(l => zeros(Float64, n, n) for l in 0:5)
    factors = Dict{Int, Any}()
    
    b_buf = zeros(Float64, n)
    y_buf = zeros(Float64, n)
    y_orb_buffer = zeros(Float64, n)
    y_total_buffer = zeros(Float64, n)

    scratch_source = zeros(Float64, n)
    scratch_y = zeros(Float64, n)
    rk_cache = Dict{UInt64, Float64}()

    return SolverWorkspace{K}(
        basis, S, T, V, R, R2, R_inv2, R_inv3, J, 
        K_mats, F_mats, tensors, factors, 
        zeros(Float64, K), zeros(Float64, K),
        b_buf, y_buf, y_orb_buffer, y_total_buffer, scratch_source, scratch_y, rk_cache
    )
end


function assemble_geometry(basis::BSplineBasis{K}, Z::Float64; calc_R_matrices::Bool=true) where {K}
    n = basis.num_splines

    S  = zeros(Float64, n, n)
    T  = zeros(Float64, n, n)
    V  = zeros(Float64, n, n)
    R  = zeros(Float64, n, n)   
    
    R2 = calc_R_matrices ? zeros(Float64, n, n) : zeros(Float64, 0, 0)
    R_inv2 = calc_R_matrices ? zeros(Float64, n, n) : zeros(Float64, 0, 0)
    R_inv3 = calc_R_matrices ? zeros(Float64, n, n) : zeros(Float64, 0, 0)

    tensors = Vector{Array{Float64, 3}}(undef, length(basis.knots) - 1)
    
    vals = zeros(Float64, K)
    derivs = zeros(Float64, K)

    for i in 1:(length(basis.knots) - 1)
        t_a = basis.knots[i]; t_b = basis.knots[i+1]
        if t_a == t_b; continue; end
        
        mid = (t_a + t_b)/2; scale = (t_b - t_a)/2
        first_global = i - K + 1

        W_local = zeros(Float64, K, K, K)

        for q in 1:length(basis.gl_nodes)
            r = scale * basis.gl_nodes[q] + mid
            r2 = r * r
            w = scale * basis.gl_weights[q]
            
            inv_r = (r > 1e-12) ? 1.0/r : 0.0
            inv_r2 = inv_r * inv_r
            inv_r3 = inv_r * inv_r2

            eval_bspline_kernel!(vals, derivs, Val(true), Val(true), i, r, basis.knots, Val(K))

            for a in 1:K
                g_a = first_global + a - 1
                Na = vals[a]; dNa = derivs[a]
                wa_tensor_factor = Na * w * inv_r 

                for b in a:K
                    g_b = first_global + b - 1
                    Nb = vals[b]; dNb = derivs[b]
                    
                    if g_a >= 1 && g_a <= n && g_b >= 1 && g_b <= n
                        term_S = w * Na * Nb
                        
                        S[g_a, g_b]  += term_S
                        T[g_a, g_b]  += w * 0.5 * dNa * dNb
                        V[g_a, g_b]  -= Z * term_S * inv_r
                        R_inv2[g_a, g_b] += term_S * inv_r2
                        
                        if calc_R_matrices
                            R[g_a, g_b]      += term_S * r
                            R2[g_a, g_b]     += term_S * r2
                            R_inv3[g_a, g_b] += term_S * inv_r3
                        end
                    end
                    
                    term_tensor_base = wa_tensor_factor * Nb
                    for k in 1:K
                        val = vals[k] * term_tensor_base
                        W_local[k, a, b] += val
                        if a != b
                            W_local[k, b, a] += val
                        end
                    end
                end
            end
        end
        tensors[i] = W_local
    end
    
    if calc_R_matrices
        return Symmetric(S), Symmetric(T), Symmetric(V), Symmetric(R), Symmetric(R2), Symmetric(R_inv2), Symmetric(R_inv3), tensors
    else
        return Symmetric(S), Symmetric(T), Symmetric(V), Symmetric(R), R2, R_inv2, R_inv3, tensors
    end
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


# Overload assemble_J_matrix! to accept a destination matrix
function assemble_J_matrix!(ws::SolverWorkspace{K}, dest_mat::Matrix{Float64}, y_coeffs) where {K}
    n = ws.basis.num_splines
    fill!(dest_mat, 0.0)

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
                    dest_mat[g_a, g_b] += val_J
                end
            end
        end
    end
    
    return dest_mat
end

# Build a specific J matrix for a single orbital and a specific scaling factor
function build_specific_J_matrix!(ws::SolverWorkspace, dest_mat::Matrix{Float64}, orb::Orbital, scale::Float64)
    fill!(ws.y_orb_buffer, 0.0)
    
    # Solve Poisson for the charge density of the provided orbital (monopole k=0 is implicit for J)
    solve_poisson_J!(ws, ws.y_orb_buffer, orb)
    
    # Scale the potential by the provided specific factor (e.g., w-1)
    ws.y_total_buffer .= scale .* ws.y_orb_buffer
    
    return assemble_J_matrix!(ws, dest_mat, ws.y_total_buffer)
end

# Build a specific K matrix bypassing the multipole logic rules
function build_specific_K_matrix!(ws::SolverWorkspace{K}, dest_mat::Matrix{Float64}, orb::Orbital, target_k::Int, coeff::Float64) where {K}
    fill!(dest_mat, 0.0)
    n = ws.basis.num_splines
    
    b_vec = zeros(Float64, n)
    y_k   = zeros(Float64, n)
    
    for nu in 1:n
        fill!(b_vec, 0.0)
        
        # --- Build Source Vector ---
        for i in 1:(length(ws.basis.knots)-1)
            if !isassigned(ws.interaction_tensors, i); continue; end
            W = ws.interaction_tensors[i]
            first = i - K + 1
            
            loc_nu = nu - first + 1
            if loc_nu < 1 || loc_nu > K; continue; end
            
            c_loc = MVector{K, Float64}(undef)
            for idx in 1:K
                g = first + idx - 1
                c_loc[idx] = (g >= 1 && g <= n) ? orb.coeffs[g] : 0.0
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
        
        # --- Solve for the targeted k multipole ---
        solve_generalized_poisson!(ws, y_k, b_vec, target_k)
        
        # --- Accumulate with the specific coefficient ---
        for i in 1:(length(ws.basis.knots)-1)
            if !isassigned(ws.interaction_tensors, i); continue; end
            W = ws.interaction_tensors[i]
            first = i - K + 1

            c_psi = MVector{K, Float64}(undef)
            c_y   = MVector{K, Float64}(undef)
            
            for idx in 1:K
                g = first + idx - 1
                if g >= 1 && g <= n
                    c_psi[idx] = orb.coeffs[g]
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
                dest_mat[g_mu, nu] += val * coeff
            end
        end
    end 
    
    # Symmetrize the specific K matrix
    for i in 1:n
        for j in (i+1):n
            dest_mat[j, i] = dest_mat[i, j]
        end
    end
    
    return Symmetric(dest_mat)
end
