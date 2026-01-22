using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using LinearAlgebra
using Plots
using Printf

function solve_helium()
    println("\n================================================================")
    println("   EXPERIMENTO: HELIO Y EL EFECTO DE APANTALLAMIENTO (V2)")
    println("================================================================\n")

    # 1. Definir el sistema
    helium = Atom(2.0, [Orbital(1, 0, 2.0)]) 
    orb1s = helium.orbitals[1]

    # 2. Configuración de la Base
    R_MAX = 1.10644808 # Critical radius
    N_ELEMS = 80      
    ORDER = 5
    GAMMA = 2.0
    basis = generate_basis(R_MAX, N_ELEMS, ORDER, γ=GAMMA)
    
    println("A. > Base generada con $(basis.num_splines) splines.")

    # 3. Operadores Base
    T, V_nuc, S = assemble_core(basis, helium.Z)

    # --------------------------------------------------------------------------
    # B. GUESS INICIAL (He+ sin interacción)
    # --------------------------------------------------------------------------
    println("\nB. Calculando Estado Inicial (Sin interacción)...")

    # Usamos la versión simple de solve_orbital!
    solve_orbital!(orb1s, helium, basis)

    # Guardamos copia para la gráfica final (Estado hidrogenoide)
    coeffs_initial = copy(orb1s.coeffs)
    E_initial_total = orb1s.energy * 2.0
    
    println("   > Energía Inicial (He+ x 2): $(E_initial_total) Ha")   

    # --------------------------------------------------------------------------
    # C. CICLO SCF (Self-Consistent Field)
    # --------------------------------------------------------------------------
    println("\nC. Iniciando SCF (Encendiendo repulsión e-e)...")
    
    MAX_ITER = 30
    MIXING = 0.6 # Factor de mezcla
    history_E = Float64[] 
    E_old = E_initial_total

    for iter in 1:MAX_ITER
        # 1. Potencial de Hartree (J)
        y_coeffs = solve_poisson_potential(basis, orb1s.coeffs, T)
        J_mat = assemble_J_matrix(basis, y_coeffs)

        # 2. Guardar coeficientes previos para el mixing
        old_coeffs = copy(orb1s.coeffs)

        # 3. Resolver usando la versión SCF de solve_orbital! (pasando J_mat)
        solve_orbital!(orb1s, helium, basis, J_mat)

        # 4. Energía Total Correcta: E = 2*eps - <J>
        E_J = dot(orb1s.coeffs, J_mat * orb1s.coeffs)
        E_total = 2 * orb1s.energy - E_J

        # 5. Mixing para mejorar convergencia
        orb1s.coeffs = MIXING * orb1s.coeffs + (1.0 - MIXING) * old_coeffs

        # 6. Re-normalizar debido al mixing
        n = sqrt(dot(orb1s.coeffs, S * orb1s.coeffs))
        orb1s.coeffs ./= n

        push!(history_E, E_total)
        
        diff = abs(E_total - E_old)
        @printf("   Iter %2d: E = %.8f Ha (Δ = %.1e)\n", iter, E_total, diff)
        
        if diff < 1e-8
            println("   >> Convergencia alcanzada.")
            break
        end
        
        E_old = E_total
    end

    # --------------------------------------------------------------------------
    # D. ANÁLISIS VISUAL
    # --------------------------------------------------------------------------
    println("\nD. Generando Gráficas...")
    
    r_plot = range(0, 5.0, length=200)
    
    # Usamos eval_expansion con los coeffs guardados
    rho_initial = [eval_expansion(coeffs_initial, basis, r)^2 for r in r_plot]
    rho_final   = [eval_expansion(orb1s.coeffs, basis, r)^2 for r in r_plot]

    # Gráfica 1: Energía
    p1 = plot(history_E, marker=:circle, label="Energía SCF",
              xlabel="Iteración", ylabel="E (Ha)", title="Convergencia SCF")
    hline!(p1, [-2.86168], label="Límite HF", linestyle=:dash, color=:red)

    # Gráfica 2: Apantallamiento
    p2 = plot(r_plot, rho_initial, label="He+ (Hidrogenoide)", color=:blue, ls=:dash)
    plot!(p2, r_plot, rho_final, label="He (SCF)", color=:red, lw=2, fill=(0, 0.1, :red))
    
    plot!(p2, title="Apantallamiento Electrónico", xlabel="r (u.a.)", ylabel="|P(r)|²")

    l = @layout [a; b]
    p_final = plot(p1, p2, layout=l, size=(700, 800))
    display(p_final)
    savefig(p_final, "helio_apantallamiento.png")
end

solve_helium()
