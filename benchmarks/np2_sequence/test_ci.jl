using Pkg
Pkg.activate("../..")
using WignerSymbols

function rme_C(l1, k, l2)
    return (-1)^l1 * sqrt((2*l1+1)*(2*l2+1)) * wigner3j(Float64, l1, k, l2, 0, 0, 0)
end

function off_diag_l2_l2(l1, l2, L, S, k)
    w6j = wigner6j(Float64, l1, l1, L, l2, l2, k)
    return (-1)^L * w6j * rme_C(l1, k, l2)^2
end

println("3P (L=1, S=1):")
println("k=1: ", off_diag_l2_l2(1, 2, 1, 1, 1))
println("k=3: ", off_diag_l2_l2(1, 2, 1, 1, 3))

println("1D (L=2, S=0):")
println("k=1: ", off_diag_l2_l2(1, 2, 2, 0, 1))
println("k=3: ", off_diag_l2_l2(1, 2, 2, 0, 3))
