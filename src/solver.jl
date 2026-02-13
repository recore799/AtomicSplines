using LinearAlgebra
using SparseArrays
using Arpack

"""
    solve_eigen(H, S; nev=1, method=:auto)

Solves the generalized eigenvalue problem Hc = ESc.

Arguments:
- `H`, `S`: AbstractMatrices (Dense or Sparse).
- `nev`: Number of eigenvalues requested.
- `method`: 
    - `:dense` forces use of LAPACK `eigen` (returns all states).
    - `:arnoldi` forces use of ARPACK `eigs` (returns `nev` states).
    - `:auto` chooses based on matrix size (cutoff N=500).
"""
function solve_eigen(H::AbstractMatrix, S::AbstractMatrix; nev=1, method=:auto)
    N = size(H, 1)
    
    # 1. Handle Boundary Conditions (Dirichlet)
    # We slice the matrices to remove the first and last rows/cols
    inner = 2:(N-1)
    H_in = H[inner, inner]
    S_in = S[inner, inner]
    
    # Determine method
    if method == :auto
        use_dense = (size(H_in, 1) < 500)
    else
        use_dense = (method == :dense)
    end

    # 2. SOLVE
    if use_dense
        # --- DENSE SOLVER (LAPACK) ---
        vals, vecs = eigen(Matrix(H_in), Matrix(S_in))
        
        # Select only the ones we asked for
        sel_vals = vals[1:nev]
        sel_vecs = vecs[:, 1:nev]
        
    else
        # --- SPARSE SOLVER (ARPACK) ---
        # which=:SR means Smallest Real
        try
            sel_vals, sel_vecs = eigs(H_in, S_in; nev=nev, which=:SR)
            
            # Arpack sometimes returns Complex numbers with 0.0im imaginary part.
            sel_vals = real.(sel_vals)
            sel_vecs = real.(sel_vecs)
        catch e
            println("Arpack failed to converge. Falling back to dense solver.")
            return solve_eigen(H, S; nev=nev, method=:dense)
        end
    end

    # 3. RECONSTRUCT (Add Zeros to boundaries)
    full_vecs = zeros(N, nev)
    full_vecs[inner, :] = sel_vecs
    
    return sel_vals, full_vecs
end


