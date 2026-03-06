using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using JLD2
using LinearAlgebra
using Printf
using WignerSymbols

struct Config1S
    orb1::Orbital
    orb2::Orbital
end

function calc_h1_1S(ws, a::Orbital, b::Orbital, c::Orbital, d::Orbital)
    # The 1-electron operator requires identical angular momentum
    if a.l != c.l
        return 0.0
    end
    
    # B-spline eigenvectors are S-orthogonal. Calculate exact numerical overlaps.
    # Since config pairs always have orb1.l == orb2.l, a, b, c, d all share the same l here.
    S_bd = dot(b.coeffs, ws.S * d.coeffs)
    S_ac = dot(a.coeffs, ws.S * c.coeffs)
    S_bc = dot(b.coeffs, ws.S * c.coeffs)
    S_ad = dot(a.coeffs, ws.S * d.coeffs)
    
    # Self-overlaps to determine normalization
    S_ab = dot(a.coeffs, ws.S * b.coeffs)
    S_cd = dot(c.coeffs, ws.S * d.coeffs)
    
    N_ab = (a.l == b.l && abs(S_ab) > 0.99) ? sqrt(2.0) : 1.0
    N_cd = (c.l == d.l && abs(S_cd) > 0.99) ? sqrt(2.0) : 1.0
    
    # Add the Centrifugal Barrier
    H_core_a = ws.T .+ ws.V .+ ((a.l * (a.l + 1)) / 2.0) .* ws.V2
    H_core_b = ws.T .+ ws.V .+ ((b.l * (b.l + 1)) / 2.0) .* ws.V2
    
    h_ac = dot(a.coeffs, H_core_a * c.coeffs)
    h_bd = dot(b.coeffs, H_core_b * d.coeffs)
    h_ad = dot(a.coeffs, H_core_a * d.coeffs)
    h_bc = dot(b.coeffs, H_core_b * c.coeffs)
    
    term1 = h_ac * S_bd + h_bd * S_ac
    term2 = h_ad * S_bc + h_bc * S_ad
    
    return (term1 + term2) / (N_ab * N_cd)
end

function calc_h2_1S(ws, a::Orbital, b::Orbital, c::Orbital, d::Orbital)
    l = a.l
    l_prime = c.l
    
    # Fix normalization here as well
    S_ab = dot(a.coeffs, ws.S * b.coeffs)
    S_cd = dot(c.coeffs, ws.S * d.coeffs)
    
    N_ab = (a.l == b.l && abs(S_ab) > 0.99) ? sqrt(2.0) : 1.0
    N_cd = (c.l == d.l && abs(S_cd) > 0.99) ? sqrt(2.0) : 1.0
    
    val = 0.0
    k_min = abs(l - l_prime)
    k_max = l + l_prime
    
    for k in k_min:2:k_max
        w3j = wigner3j(Float64, l, k, l_prime, 0, 0, 0)
        
        if abs(w3j) > 1e-10
            C_k = ((-1)^k) * sqrt((2l + 1) * (2l_prime + 1)) * (w3j^2)
            
            rk1 = compute_Rk(ws, a, b, c, d, k)
            rk2 = compute_Rk(ws, a, b, d, c, k)
            
            val += C_k * (rk1 + rk2)
        end
    end
    
    return val / (N_ab * N_cd)
end

function build_helium_ci_matrix(ws, orbitals, virtuals)
    orb_1s = orbitals[1]
    ci_basis = Config1S[]
    
    # 1. Ground State (1s^2)
    push!(ci_basis, Config1S(orb_1s, orb_1s))
    
    # 2. Singly Excited States (1s ns)
    for v in virtuals
        if v.l == 0
            push!(ci_basis, Config1S(orb_1s, v))
        end
    end
    
    # 3. Doubly Excited States (n1_l n2_l)
    for i in 1:length(virtuals)
        for j in i:length(virtuals)
            # Both electrons must jump into the same angular momentum shell to maintain L=0
            if virtuals[i].l == virtuals[j].l
                push!(ci_basis, Config1S(virtuals[i], virtuals[j]))
            end
        end
    end
    
    N = length(ci_basis)
    H_CI = zeros(Float64, N, N)
    println("Construyendo matriz CI COMPLETA para Helio ($N x $N)...")
    
    for i in 1:N
        for j in i:N
            cI = ci_basis[i]
            cJ = ci_basis[j]
            
            h1 = calc_h1_1S(ws, cI.orb1, cI.orb2, cJ.orb1, cJ.orb2)
            h2 = calc_h2_1S(ws, cI.orb1, cI.orb2, cJ.orb1, cJ.orb2)
            
            val = h1 + h2
            H_CI[i, j] = val
            H_CI[j, i] = val
        end
    end
    
    return Symmetric(H_CI)
end

function run_helium_ci(filename)
    println("=== Iniciando Configuration Interaction ===")
    println("Cargando datos de $filename...")
    
    # 1. Load the converged data
    @load filename orbitals all_virtuals E_total R_max
    
    # 2. Reconstruct the computational workspace
    println("Reconstruyendo el espacio de trabajo B-spline...")
    N_elems = 200 # Must match exactly what you used in solve_helium!
    Z = 2.0
    basis = generate_basis(R_max, N_elems, Val(7), γ=2.5)
    ws = init_scf_workspace(basis, Z)
    
    # 3. Build the Hamiltonian
    H_CI = build_helium_ci_matrix(ws, orbitals, all_virtuals)
    
    # 4. Diagonalize
    println("Diagonalizando la matriz CI...")
    t0 = time()
    evals, evecs = eigen(H_CI)
    elapsed = time() - t0
    
    # 5. Extract the physics
    E_HF = E_total
    E_CI = evals[1]  # The new, correlated ground state
    E_corr = E_CI - E_HF
    
    println("\n=== Resultados Finales ===")
    @printf("Tiempo de diag.      : %.4f s\n", elapsed)
    @printf("Energía Hartree-Fock : %12.6f Ha\n", E_HF)
    @printf("Energía CI (Total)   : %12.6f Ha\n", E_CI)
    println("-"^30)
    @printf("Energía Correlación  : %12.6f Ha\n", E_corr)
    println("==========================\n")
    
    return evals, evecs, H_CI
end

# Run it!
run_helium_ci("helium_results_R20.0.jld2")
