using Pkg
Pkg.activate("../..")
using WignerSymbols

println("w6j(1,1,1,2,2,1) = ", wigner6j(Float64, 1, 1, 1, 2, 2, 1))
println("w6j(1,1,1,2,2,3) = ", wigner6j(Float64, 1, 1, 1, 2, 2, 3))
println("w6j(1,1,2,2,2,1) = ", wigner6j(Float64, 1, 1, 2, 2, 2, 1))
println("w6j(1,1,2,2,2,3) = ", wigner6j(Float64, 1, 1, 2, 2, 2, 3))
