using AtomicSplines
using Test
using LinearAlgebra

function solve_hydrogen()
    # Grid: Large box, exponential spacing
    R_MAX = 60.0
    N_ELEMS = 100
    ORDER = 6
    basis = generate_basis(R_MAX, N_ELEMS, ORDER, γ=3.0)
    
    # Z=1 for Hydrogen
    # This line does the heavy lifting from the package
    T, V, S = assemble_core(basis, 1.0)
    
    # Hamiltonian: H = T + V
    H = T + V
    
    # Boundary Conditions (Dirichlet at r=0 and r=R_max)
    # We strip the first and last B-splines.
    inner = 2:(basis.num_splines - 1)
    
    # Generalized Eigenvalue Problem: H c = E S c
    # We must convert sparse to dense for the 'eigen' solver in tests
    # (Arnoldi methods are better for production, but this is robust for testing)
    H_inner = Matrix(H[inner, inner])
    S_inner = Matrix(S[inner, inner])
    
    vals, _ = eigen(H_inner, S_inner)
    
    return vals[1] # The Ground State Energy
end

@testset "AtomicSplines.jl Tests" begin

    @testset "1. Basis Properties" begin
        # Check that we have the right number of functions
        # N_elements + Order - 1 (Standard B-Spline count)
        b = generate_basis(10.0, 10, 4)
        @test b.num_splines == 13
        @test b.knots[1] == 0.0
        @test b.knots[end] == 10.0
    end

    @testset "2. Physics: Hydrogen Atom" begin
        # Exact Energy: -0.5 Ha
        
        energy = solve_hydrogen()
        
        println("\n   > Hydrogen Ground State: $energy Ha")
        println("   > Exact Target:        -0.5 Ha")
        
        @test isapprox(energy, -0.5, atol=1e-5)
    end
end
