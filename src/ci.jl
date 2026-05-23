struct Term
    L::Int
    S::Float64
end

function Base.show(io::IO, t::Term)
    multiplicity = Int(2 * t.S + 1)
    L_chars = ['S', 'P', 'D', 'F', 'G', 'H', 'I']
    L_str = t.L < length(L_chars) ? L_chars[t.L + 1] : string(t.L)
    print(io, "^$(multiplicity)$(L_str)")
end

@inline function get_orb_id(o::Orbital)::UInt64
    return (UInt64(o.n) << 4) | UInt64(o.l)
end

@inline function pack_rk_key(a::Orbital, b::Orbital, c::Orbital, d::Orbital, k::Int)::UInt64
    id_a = get_orb_id(a)
    id_c = get_orb_id(c)
    id_b = get_orb_id(b)
    id_d = get_orb_id(d)

    id_a, id_c = id_a < id_c ? (id_a, id_c) : (id_c, id_a)
    id_b, id_d = id_b < id_d ? (id_b, id_d) : (id_d, id_b)

    pair1 = (id_a << 14) | id_c
    pair2 = (id_b << 14) | id_d

    pair1, pair2 = pair1 < pair2 ? (pair1, pair2) : (pair2, pair1)

    return (UInt64(k) << 56) | (pair1 << 28) | pair2
end

function get_cached_Rk!(ws::SolverWorkspace, a::Orbital, b::Orbital, c::Orbital, d::Orbital, k::Int)
    key = pack_rk_key(a, b, c, d, k) 
    
    if haskey(ws.rk_cache, key)
        return ws.rk_cache[key]
    end
    
    val = compute_Rk(ws, a, b, c, d, k)
    ws.rk_cache[key] = val
    return val
end

const CFP_P_SHELL = Dict{Int, Dict{Term, Vector{Tuple{Term, Float64}}}}(
    2 => Dict(
        Term(1, 1.0) => [(Term(1, 0.5), 1.0)],          
        Term(2, 0.0) => [(Term(1, 0.5), 1.0)],          
        Term(0, 0.0) => [(Term(1, 0.5), 1.0)]           
    ),
    3 => Dict(
        Term(0, 1.5) => [ 
            (Term(1, 1.0), 1.0)                         
        ],
        Term(2, 0.5) => [ 
            (Term(1, 1.0), -sqrt(1.0/2.0)),             
            (Term(2, 0.0), sqrt(1.0/2.0))               
        ],
        Term(1, 0.5) => [ 
            (Term(1, 1.0), -sqrt(1.0/2.0)),          
            (Term(0, 0.0), sqrt(5.0/18.0)),             
            (Term(2, 0.0), -sqrt(4.0/9.0))              
        ]
    )
)

function extract_virtuals(evals::Vector{Float64}, evecs::Matrix{Float64}, l::Int, num_occupied::Int, N_virt::Int, active_idx::UnitRange{Int}, n_splines::Int, ws; offset::Int=0)
    virtuals = Orbital[]
    start_idx = num_occupied + 1 + offset
    end_idx = min(start_idx + N_virt - 1, length(evals))
    
    for i in start_idx:end_idx
        pseudo_n = l + 1 + num_occupied + (i - start_idx) 
        
        virt_orb = Orbital(pseudo_n, l, 0.0)
        virt_orb.energy = evals[i]
        
        virt_orb.coeffs = zeros(Float64, n_splines)
        virt_orb.coeffs[active_idx] = evecs[:, i]
        
        norm_factor = sqrt(dot(virt_orb.coeffs, ws.S * virt_orb.coeffs))
        virt_orb.coeffs ./= norm_factor
        
        push!(virtuals, virt_orb)
    end
    
    return virtuals
end

function get_cfp_expansion(l::Int, k::Int, target_term::Term)
    if k == 1
        return [(Term(0, 0.0), 1.0)]
    end
    
    if l == 1
        if haskey(CFP_P_SHELL, k) && haskey(CFP_P_SHELL[k], target_term)
            return CFP_P_SHELL[k][target_term]
        else
            error("State $(target_term) is forbidden by Pauli principle for p^$(k).")
        end
    else
        error("Only p-shell CFPs are implemented right now.")
    end
end

function evaluate_equivalent_shell_energy(l::Int, k::Int, target::Term)
    parents = get_cfp_expansion(l, k, target)
    total_energy = 0.0
    
    for (parent_term, cfp_weight) in parents
        interaction_energy = compute_wigner_interaction(parent_term, target, l)
        total_energy += (cfp_weight^2) * interaction_energy
    end
    
    num_pairs = k * (k - 1) / 2.0
    return num_pairs * total_energy
end

function compute_Rk(ws::SolverWorkspace{K}, a::Orbital, b::Orbital, c::Orbital, d::Orbital, k_mult::Int) where {K}
    n = ws.basis.num_splines
    
    source = ws.scratch_source
    y_k    = ws.scratch_y
    fill!(source, 0.0)
    fill!(y_k, 0.0)

    @inbounds for i in 1:(length(ws.basis.knots)-1)
        if !isassigned(ws.interaction_tensors, i); continue; end
        W = ws.interaction_tensors[i]
        first = i - K + 1

        c_b = MVector{K, Float64}(undef)
        c_d = MVector{K, Float64}(undef)
        for idx in 1:K
            g = first + idx - 1
            valid = (1 <= g <= n)
            c_b[idx] = valid ? b.coeffs[g] : 0.0
            c_d[idx] = valid ? d.coeffs[g] : 0.0
        end
        
        src_local = zeros(MVector{K, Float64})
        
        for beta in 1:K
            for alpha in 1:K
                term = c_b[alpha] * c_d[beta]
                @simd for k_idx in 1:K
                    src_local[k_idx] += term * W[k_idx, alpha, beta]
                end
            end
        end
       
        for k_idx in 1:K
            g_k = first + k_idx - 1
            if 1 <= g_k <= n
                source[g_k] += src_local[k_idx]
            end
        end
    end
    
    solve_generalized_poisson!(ws, y_k, source, k_mult)
    
    Rk_value = 0.0
    
    @inbounds for i in 1:(length(ws.basis.knots)-1)
        if !isassigned(ws.interaction_tensors, i); continue; end
        W = ws.interaction_tensors[i]
        first = i - K + 1

        c_a = MVector{K, Float64}(undef)
        c_c = MVector{K, Float64}(undef)
        c_y = MVector{K, Float64}(undef)
        
        for idx in 1:K
            g = first + idx - 1
            valid = (1 <= g <= n)
            c_a[idx] = valid ? a.coeffs[g] : 0.0
            c_c[idx] = valid ? c.coeffs[g] : 0.0
            c_y[idx] = valid ? y_k[g] : 0.0
        end

        val_local = 0.0
        
        for gamma in 1:K
            for beta in 1:K
                term = c_c[beta] * c_y[gamma]
                @simd for alpha in 1:K
                    val_local += c_a[alpha] * term * W[alpha, beta, gamma]
                end
            end
        end
        Rk_value += val_local
    end
    
    return Rk_value
end
