using Pkg
Pkg.activate(".")
include("benchmarks/np2_sequence/tin_rohf.jl")

function test_grid(N, g)
    println("Testing N_elems = $N, gamma = $g")
    
    # We will basically copy the solve_tin_rohf body but just do 10 iterations to see the energy trajectory,
    # or just copy the file and replace N_elems and gamma.
end
