using LinearAlgebra
using Printf
using FastGaussQuadrature # Faster than standard GaussLegendre if available, else use standard

# ==============================================================================
# 1. CORE B-SPLINE ENGINE (Optimized)
# ==============================================================================

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

# ==============================================================================
# 2. MATRIX ASSEMBLERS (Element-Wise / Vectorized)
# ==============================================================================

"""
Builds T (Kinetic), V_nuc (Nuclear), and S (Overlap) in a single pass.
Iterates over ELEMENTS (intervals), not functions.
"""
function assemble_core(basis, Z)
    n = basis.num_splines
    k = basis.order
    
    # Use dense matrices for N < 500. 
    T = zeros(Float64, n, n)
    V = zeros(Float64, n, n)
    S = zeros(Float64, n, n)
    
    scratch = BSplineScratch(Vector{Float64}(undef, k), Vector{Float64}(undef, k))
    
    # Gauss-Legendre Quadrature
    gl_p, gl_w = gausslegendre(k + 2) 

    # Iterate over knot intervals (elements)
    # The active range of splines is roughly [1, n+k] in the knot vector
    for i in 1:(length(basis.knots) - 1)
        t_a = basis.knots[i]; t_b = basis.knots[i+1]
        if t_a == t_b; continue; end

        # Map bounds
        mid = (t_a + t_b) / 2
        scale = (t_b - t_a) / 2
        
        # Identify active B-splines: B_{i-k+1} ... B_{i} are non-zero here
        first_global = i - k + 1

        for q in 1:length(gl_p)
            r = scale * gl_p[q] + mid
            w = scale * gl_w[q]

            # Kernel Call: Get Vals AND Derivs
            eval_bspline_kernel!(scratch.vals, scratch.derivs, Val(true), Val(true), 
                                 i, k, r, basis.knots)

            # Add to matrix blocks
            for a in 1:k
                g_a = first_global + a - 1
                if g_a < 1 || g_a > n; continue; end
                
                # Pre-fetch for inner loop
                Na  = scratch.vals[a]
                dNa = scratch.derivs[a]

                for b in 1:k
                    g_b = first_global + b - 1
                    if g_b < 1 || g_b > n; continue; end
                    
                    Nb  = scratch.vals[b]
                    dNb = scratch.derivs[b]
                    
                    # Overlap: <a|b>
                    S[g_a, g_b] += w * Na * Nb
                    
                    # Kinetic: 0.5 * <da|db>
                    T[g_a, g_b] += w * 0.5 * dNa * dNb
                    
                    # Nuclear: -Z * <a|1/r|b>
                    # Avoid singularity at r=0 (though GL points shouldn't hit 0)
                    if r > 1e-12
                        V[g_a, g_b] += w * (-Z / r) * Na * Nb
                    end
                end
            end
        end
    end
    
    return Symmetric(T), Symmetric(V), Symmetric(S)
end

"""
Builds the Hartree Potential Matrix J.
Calculates density y(r) on the fly without extra searches.
"""
function assemble_J_matrix(basis, y_coeffs)
    n = basis.num_splines; k = basis.order
    J_mat = zeros(Float64, n, n)
    
    scratch = BSplineScratch(Vector{Float64}(undef, k), Vector{Float64}(undef, k))
    gl_p, gl_w = gausslegendre(k + 2)

    for i in 1:(length(basis.knots) - 1)
        t_a = basis.knots[i]; t_b = basis.knots[i+1]
        if t_a == t_b; continue; end
        
        first_global = i - k + 1
        mid = (t_a + t_b)/2; scale = (t_b - t_a)/2

        for q in 1:length(gl_p)
            r = scale * gl_p[q] + mid
            w = scale * gl_w[q]

            # Need Values only
            eval_bspline_kernel!(scratch.vals, scratch.derivs, Val(true), Val(false), 
                                 i, k, r, basis.knots)

            # Compute y(r) = r * V_H(r)
            # Actually, usually y_coeffs ARE the density projection coeffs.
            # Let's assume y_coeffs represents the solution to Poisson: U(r) = y(r)/r
            # So here we reconstruct y(r)
            y_val = 0.0
            for idx in 1:k
                g_idx = first_global + idx - 1
                if g_idx >= 1 && g_idx <= n
                    y_val += y_coeffs[g_idx] * scratch.vals[idx]
                end
            end
            
            # Potential V_H = y(r) / r
            if r < 1e-12; continue; end # Avoid singularity
            V_H = y_val / r
            
            # Fill Matrix J_ab = <a | V_H | b>
            for a in 1:k
                g_a = first_global + a - 1
                if g_a < 1 || g_a > n; continue; end
                
                term_a = scratch.vals[a] * w * V_H
                
                for b in 1:k
                    g_b = first_global + b - 1
                    if g_b < 1 || g_b > n; continue; end
                    
                    J_mat[g_a, g_b] += term_a * scratch.vals[b]
                end
            end
        end
    end
    return Symmetric(J_mat)
end

function solve_poisson_potential(basis, orbital_coeffs, T_matrix)
    n = basis.num_splines
    k = basis.order
    
    # The RHS vector b
    b = zeros(Float64, n)
    
    scratch = BSplineScratch(Vector{Float64}(undef, k), Vector{Float64}(undef, k))
    # Use slightly higher quadrature for the Coulomb term (degree ~3k)
    gl_p, gl_w = gausslegendre(k + 4)

    for i in 1:(length(basis.knots) - 1)
        t_a = basis.knots[i]; t_b = basis.knots[i+1]
        if t_a == t_b; continue; end
        
        first_global = i - k + 1
        mid = (t_a + t_b)/2; scale = (t_b - t_a)/2

        for q in 1:length(gl_p)
            r = scale * gl_p[q] + mid
            w = scale * gl_w[q]
            
            # 1. Eval basis
            eval_bspline_kernel!(scratch.vals, scratch.derivs, Val(true), Val(false), 
                                 i, k, r, basis.knots)
            
            # 2. Reconstruct radial function u(r)
            u_val = 0.0
            for idx in 1:k
                g_idx = first_global + idx - 1
                if g_idx >= 1 && g_idx <= n
                    u_val += orbital_coeffs[g_idx] * scratch.vals[idx]
                end
            end
            
            # 3. Source Term: -y'' = u^2 / r
            #    We solve K*y = b. K is positive definite (integral of B' B').
            #    The operator is -d^2/dr^2. So LHS is positive.
            #    RHS must be positive: u^2/r.
            if r > 1e-12
                source = u_val^2 / r
            else
                source = 0.0 
            end
            
            # 4. Fill b vector
            for a in 1:k
                g_a = first_global + a - 1
                if g_a < 1 || g_a > n; continue; end
                b[g_a] += w * scratch.vals[a] * source
            end
        end
    end
    
    # 5. Solve K * y = b
    #    K = integral(B' B') = 2 * T_matrix (stiffness)
    K = 2.0 .* T_matrix
    
    # --- BOUNDARY CONDITIONS ---
    # y(0) = 0 (Potential is finite, y=rV -> 0)
    # y(R) = 1 (Total Charge Enclosed = 1 electron)
    #
    # We solve for inner nodes [2:n-1].
    # The equation for node (n-1) involves y_n:
    # K[n-1, ...] * y_... + K[n-1, n] * y_n = b[n-1]
    # Move known boundary term to RHS:
    
    inner = 2:(n-1)
    y_boundary = 1.0 
    
    # Adjust RHS
    b_inner = b[inner] - K[inner, n] * y_boundary
    
    # Solve
    y_inner = K[inner, inner] \ b_inner
    
    # Reconstruct full vector
    y_full = zeros(Float64, n)
    y_full[inner] = y_inner
    y_full[n] = y_boundary # Enforce y(R) = 1
    
    return y_full
end

# ==============================================================================
# 4. MAIN SOLVER
# ==============================================================================

function solve_helium_optimized()
    println("================================================================")
    println("   HELIUM SOLVER (Optimized Element-Wise Engine)")
    println("================================================================")

    # 1. Setup
    R_MAX = 20.0
    N_ELEMS = 60
    ORDER = 7 # High order splines
    Z = 2.0
    
    basis = generate_basis(R_MAX, N_ELEMS, ORDER, γ=3.0)
    println("Basis: $(basis.num_splines) splines, Order $(basis.order)")

    # 2. Fixed Operators (T, V, S)
    print("Assembling Core Matrices... ")
    time_core = @elapsed T, V_nuc, S = assemble_core(basis, Z)
    println("Done ($(round(time_core, digits=4)) s)")



    # 3. Initial Guess (H_core * c = E * S * c)
    H_core = T + V_nuc
    
    # === FIX START: Apply Dirichlet BCs (u(0)=0, u(R)=0) ===
    # We slice the matrices to exclude the first and last basis functions.
    # This effectively forces c[1] = 0 and c[n] = 0.
    active = 2:(basis.num_splines - 1)
    
    # Solve generalized eigenproblem on the subspace
    evals, evecs_active = eigen(H_core[active, active], S[active, active])
    
    # Reconstruct the full vector (padding zeros)
    c_current = zeros(Float64, basis.num_splines)
    
    perm = sortperm(Real.(evals))
    c_current[active] = evecs_active[:, perm[1]]
    # === FIX END ===

    # Normalize
    norm_fac = sqrt(dot(c_current, S * c_current))
    c_current ./= norm_fac
    
    E_initial = evals[perm[1]]
    println("Initial Guess (He+): $(E_initial) Ha (Target: -2.0)")

    # 4. SCF Loop
    MAX_ITER = 30
    MIXING = 0.3
    TOL = 1e-9
    E_old = 0.0
    
    println("\nStarting SCF...")
    println("Iter | Total Energy (Ha) | Delta E    | Time (s)")
    println("--------------------------------------------------")
    
    for iter in 1:MAX_ITER
        t_start = time()

        # A. Solve Poisson for Hartree Potential
        y_coeffs = solve_poisson_potential(basis, c_current, T)
        
        # B. Build J Matrix
        J = assemble_J_matrix(basis, y_coeffs)
        
        # C. Build Fock Matrix: F = H_core + J
        #    (Factor 1.0 for J because we solve for spatial orbital occupied by 2 electrons.
        #     In RHF for He, Fock operator is h + J_1s)
        F = H_core + J
        
        # D. Solve Eigenproblem (WITH BCs)
        evals, evecs_active = eigen(F[active, active], S[active, active])

        # E. Update State
        perm = sortperm(Real.(evals))
        
        c_new = zeros(Float64, basis.num_splines)
        c_new[active] = evecs_active[:, perm[1]]
       
        # Normalize
        n_val = sqrt(dot(c_new, S * c_new))
        c_new ./= n_val
        
        # F. Calculate Total Energy
        # E_tot = 2*epsilon - <J>
        eps = evals[perm[1]]
        E_J = dot(c_new, J * c_new)
        E_total = 2 * eps - E_J
        
        # G. Mixing
        c_current = MIXING * c_current + (1.0 - MIXING) * c_new
        c_current ./= sqrt(dot(c_current, S * c_current)) # Re-norm
        
        t_iter = time() - t_start
        delta = abs(E_total - E_old)
        
        @printf("%4d | %.10f    | %.2e   | %.4f\n", iter, E_total, delta, t_iter)
        
        if delta < TOL
            println("--------------------------------------------------")
            println("Converged!")
            println("Final Energy: $(E_total) Ha")
            println("Ref Value   : -2.86168 Ha")
            break
        end
        E_old = E_total
    end
end

# Run it
solve_helium_optimized()


