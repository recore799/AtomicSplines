# Represents the L and S quantum numbers of a coupled state
struct Term
    L::Int
    S::Float64
end

# Printing format
function Base.show(io::IO, t::Term)
    multiplicity = Int(2 * t.S + 1)
    L_chars = ['S', 'P', 'D', 'F', 'G', 'H', 'I']
    L_str = t.L < length(L_chars) ? L_chars[t.L + 1] : string(t.L)
    print(io, "^$(multiplicity)$(L_str)")
end

# 1. Pack a single orbital into 14 bits
@inline function get_orb_id(o::Orbital)::UInt64
    # Ensure no type instability and pack: [ n (10 bits) | l (4 bits) ]
    return (UInt64(o.n) << 4) | UInt64(o.l)
end

# 2. Pack the whole integral into a canonical 64-bit key
@inline function pack_rk_key(a::Orbital, b::Orbital, c::Orbital, d::Orbital, k::Int)::UInt64
    id_a = get_orb_id(a)
    id_c = get_orb_id(c)
    id_b = get_orb_id(b)
    id_d = get_orb_id(d)

    # Enforce intra-electron symmetry (a < c and b < d)
    id_a, id_c = id_a < id_c ? (id_a, id_c) : (id_c, id_a)
    id_b, id_d = id_b < id_d ? (id_b, id_d) : (id_d, id_b)

    # Pack into 28-bit pairs for each electron
    pair1 = (id_a << 14) | id_c
    pair2 = (id_b << 14) | id_d

    # Enforce inter-electron symmetry (pair1 < pair2)
    pair1, pair2 = pair1 < pair2 ? (pair1, pair2) : (pair2, pair1)

    # Final 64-bit pack: [ k (8 bits) | pair1 (28 bits) | pair2 (28 bits) ]
    return (UInt64(k) << 56) | (pair1 << 28) | pair2
end

function get_cached_Rk!(ws::SolverWorkspace, a::Orbital, b::Orbital, c::Orbital, d::Orbital, k::Int)
    key = pack_rk_key(a, b, c, d, k) # (Using the UInt64 bit-packing function from earlier)
    
    if haskey(ws.rk_cache, key)
        return ws.rk_cache[key]
    end
    
    val = compute_Rk(ws, a, b, c, d, k)
    ws.rk_cache[key] = val
    return val
end


# CFP Table for the p-shell (l=1)
# Structure: k_electrons => Dict( Target_Term => Vector{Tuple{Parent_Term, Weight}} )
const CFP_P_SHELL = Dict{Int, Dict{Term, Vector{Tuple{Term, Float64}}}}(
    
    # p^2 splitting into p^1 parent
    2 => Dict(
        Term(1, 1.0) => [(Term(1, 0.5), 1.0)],          # ^3P target -> ^2P parent
        Term(2, 0.0) => [(Term(1, 0.5), 1.0)],          # ^1D target -> ^2P parent
        Term(0, 0.0) => [(Term(1, 0.5), 1.0)]           # ^1S target -> ^2P parent
    ),

    # p^3 splitting into p^2 parents
    3 => Dict(
        Term(0, 1.5) => [ # ^4S target
            (Term(1, 1.0), 1.0)                         # ^3P parent
        ],
        Term(2, 0.5) => [ # ^2D target
            (Term(1, 1.0), -sqrt(1.0/2.0)),             # ^3P parent
            (Term(2, 0.0), sqrt(1.0/2.0))               # ^1D parent
        ],
        Term(1, 0.5) => [ # ^2P target
            (Term(1, 1.0), -sqrt(1.0/2.0)),             # ^3P parent
            (Term(0, 0.0), sqrt(5.0/18.0)),             # ^1S parent
            (Term(2, 0.0), -sqrt(4.0/9.0))              # ^1D parent
        ]
    )
)

function extract_virtuals(evals::Vector{Float64}, evecs::Matrix{Float64}, l::Int, num_occupied::Int, N_virt::Int, active_idx::UnitRange{Int}, n_splines::Int, ws; offset::Int=0)
    virtuals = Orbital[]
    
    # Skip occupied states + whatever offset we need to bypass the centrifugal wall
    start_idx = num_occupied + 1 + offset
    end_idx = min(start_idx + N_virt - 1, length(evals))
    
    for i in start_idx:end_idx
        # Assign a pseudo-principal quantum number purely for bookkeeping
        pseudo_n = l + 1 + num_occupied + (i - start_idx) 
        
        # Initialize an empty virtual orbital (occ = 0.0)
        virt_orb = Orbital(pseudo_n, l, 0.0)
        virt_orb.energy = evals[i]
        
        # Map the active coefficients back to the full spline grid
        virt_orb.coeffs = zeros(Float64, n_splines)
        virt_orb.coeffs[active_idx] = evecs[:, i]
        
        # Enforce strict normalization
        norm_factor = sqrt(dot(virt_orb.coeffs, ws.S * virt_orb.coeffs))
        virt_orb.coeffs ./= norm_factor
        
        push!(virtuals, virt_orb)
    end
    
    return virtuals
end

function get_cfp_expansion(l::Int, k::Int, target_term::Term)
    if k == 1
        # Base case: a single electron has a "vacuum" parent with L=0, S=0
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
    # 1. Get the CFP expansion for our target state (e.g., p^3 ^2D)
    parents = get_cfp_expansion(l, k, target)
    
    total_energy = 0.0
    
    # 2. Sum over the parent states
    for (parent_term, cfp_weight) in parents
        
        # 3. Calculate the standard two-electron interaction 
        #    between the (k-1) parent core and the single stripped electron.
        #    (This is where you will call WignerSymbols.jl and your compute_Rk function)
        interaction_energy = compute_wigner_interaction(parent_term, target, l)
        
        # 4. Multiply by the square of the CFP weight 
        #    (squared because the bra and the ket both undergo the CFP expansion)
        total_energy += (cfp_weight^2) * interaction_energy
    end
    
    # Account for the number of pairs we can form in the shell
    num_pairs = k * (k - 1) / 2.0
    return num_pairs * total_energy
end

function compute_Rk(ws::SolverWorkspace{K}, a::Orbital, b::Orbital, c::Orbital, d::Orbital, k_mult::Int) where {K}
    n = ws.basis.num_splines
    
    source = ws.scratch_source
    y_k    = ws.scratch_y
    fill!(source, 0.0)
    fill!(y_k, 0.0)

    # --- Build Source Vector from Orbitals b and d ---
    @inbounds for i in 1:(length(ws.basis.knots)-1)
        if !isassigned(ws.interaction_tensors, i); continue; end
        W = ws.interaction_tensors[i]
        first = i - K + 1

        # Extract local coefficients for b and d
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
                # W tensor structure: W[k_idx, alpha, beta] integrates B_k * B_alpha * B_beta
                term = c_b[alpha] * c_d[beta]
                @simd for k_idx in 1:K
                    src_local[k_idx] += term * W[k_idx, alpha, beta]
                end
            end
        end
        
        # Scatter back to global source
        for k_idx in 1:K
            g_k = first + k_idx - 1
            if 1 <= g_k <= n
                source[g_k] += src_local[k_idx]
            end
        end
    end
    
    # --- Solve Generalized Poisson ---
    solve_generalized_poisson!(ws, y_k, source, k_mult)
    
    # --- Integrate a, c, and y_k ---
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

# function compute_Rk(ws::SolverWorkspace{K}, a::Orbital, b::Orbital, c::Orbital, d::Orbital, k_mult::Int) where {K}
#     n = ws.basis.num_splines
    
#     source = ws.scratch_source
#     y_k    = ws.scratch_y
#     fill!(source, 0.0)
#     fill!(y_k, 0.0)

#     # --- Build Source Vector from Orbitals b and d ---
#     for i in 1:(length(ws.basis.knots)-1)
#         if !isassigned(ws.interaction_tensors, i); continue; end
#         W = ws.interaction_tensors[i]
#         first = i - K + 1
        
#         # Extract local coefficients for b and d
#         c_b = MVector{K, Float64}(undef)
#         c_d = MVector{K, Float64}(undef)
#         for idx in 1:K
#             g = first + idx - 1
#             c_b[idx] = (g >= 1 && g <= n) ? b.coeffs[g] : 0.0
#             c_d[idx] = (g >= 1 && g <= n) ? d.coeffs[g] : 0.0
#         end
        
#         for k_idx in 1:K
#             g_k = first + k_idx - 1
#             if g_k < 1 || g_k > n; continue; end
            
#             val = 0.0
#             for alpha in 1:K
#                 # W tensor structure: W[k_idx, alpha, beta] integrates B_k * B_alpha * B_beta
#                 for beta in 1:K
#                     val += c_b[alpha] * c_d[beta] * W[k_idx, alpha, beta]
#                 end
#             end
#             source[g_k] += val
#         end
#     end
    
#     # Solve Generalized Poisson for Y^k_{bd} ---
#     # This mutates y_k directly
#     solve_generalized_poisson!(ws, y_k, source, k_mult)
    
#     # Integrate a, c, and y_k ---
#     Rk_value = 0.0
    
#     for i in 1:(length(ws.basis.knots)-1)
#         if !isassigned(ws.interaction_tensors, i); continue; end
#         W = ws.interaction_tensors[i]
#         first = i - K + 1

#         c_a = MVector{K, Float64}(undef)
#         c_c = MVector{K, Float64}(undef)
#         c_y = MVector{K, Float64}(undef)
        
#         for idx in 1:K
#             g = first + idx - 1
#             if g >= 1 && g <= n
#                 c_a[idx] = a.coeffs[g]
#                 c_c[idx] = c.coeffs[g]
#                 c_y[idx] = y_k[g]
#             else
#                 c_a[idx] = 0.0
#                 c_c[idx] = 0.0
#                 c_y[idx] = 0.0
#             end
#         end

#         for alpha in 1:K
#             for beta in 1:K
#                 for gamma in 1:K
#                     # Contract a, c, and the Y^k potential
#                     Rk_value += c_a[alpha] * c_c[beta] * c_y[gamma] * W[alpha, beta, gamma]
#                 end
#             end
#         end
#     end
    
#     return Rk_value
# end
