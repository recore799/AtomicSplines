using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using LinearAlgebra
using Printf
using Roots 

function calculate_energy(R_box, n, l)
    # 20 B-splines orden 7 (Consistente con Ting-Yun)
    n_bsplines = max(30, Int(round(3 * R_box))) 
    k_order = 7
    basis = generate_basis(R_box, n_bsplines, k_order, γ=1.0) 
    h_atom = Atom(1.0, [Orbital(n, l, 0.0)]) 
    try
        solve_orbital!(h_atom.orbitals[1], h_atom, basis)
        return h_atom.orbitals[1].energy
    catch
        return NaN
    end
end

function find_critical_radius(n, l, r_guess_min, r_guess_max)
    f(r) = calculate_energy(r, n, l)
    try
        return find_zero(f, (r_guess_min, r_guess_max))
    catch
        return NaN
    end
end

function format_val(val)
    abs_val = abs(val)
    if abs_val == 0
        return @sprintf("%.6f", val)
    elseif abs_val >= 1.0 || abs_val < 0.1
        return @sprintf("%.6e", val)
    else
        return @sprintf("%.6f", val)
    end
end

println("\n=======================================================")
println(" VALORES CALCULADOS (HIDRÓGENO CONFINADO)")
println("=======================================================\n")

# ENERGÍAS (Para comparar con Tabla 1 de Ting-Yun
radii = [10.0, 5.0, 1.0, 0.1]
states_e = [(1, 0, "1s"), (2, 1, "2p"), (3, 2, "3d")]

println("--- DATOS PARA TABLA 1 (Energías) ---")
println("r0\t\tE_1s\t\tE_2p\t\tE_3d")
println("-" ^ 60)

for r in radii
    # Imprime r0
    print(@sprintf("%.2f", r))
    
    # Calcula e imprime cada estado en la misma fila
    for (n, l, _) in states_e
        val = calculate_energy(r, n, l)
        print("\t\t" * format_val(val))
    end
    println() # Nueva línea
end
println("\n")

# RADIOS CRÍTICOS (Para comparar con Tabla 2 de Ting-Yun)
states_rc = [
    (1, 0, "1s", 1.0, 3.0),
    (2, 0, "2s", 5.0, 8.0),
    (2, 1, "2p", 4.0, 6.0),
    (3, 2, "3d", 8.0, 12.0)
]

println("--- DATOS PARA TABLA 2 (Radios Críticos) ---")
println("Estado\tr_c (Calculado)")
println("-" ^ 30)

for (n, l, label, min_r, max_r) in states_rc
    rc = find_critical_radius(n, l, min_r, max_r)
    println("$label\t" * @sprintf("%.5f", rc))
end
println("\n")

