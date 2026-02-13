using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using LinearAlgebra
using Plots
using Printf

function solve_helium_screening()
    println("\n================================================================")
    println("   EXPERIMENTO: HELIO Y EL EFECTO DE APANTALLAMIENTO (V2)")
    println("================================================================\n")

    # 1. Definir parámetros físicos
    Z = 2.0             # Carga nuclear
    R_MAX = 6.0         # Radio de caja (suficiente para el estado base)
    N_ELEMS = 60        # Elementos finitos
    ORDER = 6           # Orden de los B-splines
    GAMMA = 2.5         # Distribución de nodos
    
    # 2. Generar Base
    basis = generate_basis(R_MAX, N_ELEMS, ORDER, γ=GAMMA)
    println("A. > Base generada con $(basis.num_splines) splines.")

    # 3. Ensamblar Operadores Nucleares (T, V, S)
    print("   > Ensamblando matrices Core... ")
    T, V_nuc, S = assemble_core(basis, Z)
    println("Listo.")

    # Definir subespacio activo (Dirichlet BCs: u(0)=0, u(R)=0)
    active = 2:(basis.num_splines - 1)
    
    # --------------------------------------------------------------------------
    # B. GUESS INICIAL (He+ sin interacción)
    # --------------------------------------------------------------------------
    println("\nB. Calculando Estado Inicial (Sin interacción e-e)...")

    H_core = T + V_nuc
    
    # Resolver eigenproblema H c = E S c
    evals, evecs = eigen(H_core[active, active], S[active, active])
    
    # Reconstruir vector completo
    coeffs_initial = zeros(Float64, basis.num_splines)
    perm = sortperm(Real.(evals))
    coeffs_initial[active] = evecs[:, perm[1]]
    
    # Normalizar
    coeffs_initial ./= sqrt(dot(coeffs_initial, S * coeffs_initial))
    
    E_initial_single = evals[perm[1]]
    E_initial_total = 2 * E_initial_single # 2 electrones sin repulsión
    
    println("   > Energía Inicial (He+ x 2): $(E_initial_total) Ha")    
    
    # Inicializar ciclo SCF con el guess
    c_current = copy(coeffs_initial)

    # --------------------------------------------------------------------------
    # C. CICLO SCF (Self-Consistent Field)
    # --------------------------------------------------------------------------
    println("\nC. Iniciando SCF (Encendiendo repulsión e-e)...")
    
    MAX_ITER = 30
    MIXING = 0.3
    TOL = 1e-9
    
    history_E = Float64[] 
    E_old = E_initial_total
    
    # Guardamos energías para la gráfica
    push!(history_E, E_initial_total)

    for iter in 1:MAX_ITER
        # 1. Potencial de Hartree (J)
        #    Resuelve Poisson: V_H(r) = U(r)/r
        y_coeffs = solve_poisson_potential(basis, c_current, T)
        
        # 2. Ensamblar Matriz de Interacción Directa J
        J_mat = assemble_J_matrix(basis, y_coeffs)

        # 3. Construir Operador de Fock: F = H_core + J
        #    (Restricted Hartree-Fock para capa cerrada 1s^2)
        F = H_core + J_mat

        # 4. Resolver Eigenproblema Generalizado
        evals, evecs = eigen(F[active, active], S[active, active])
        
        # 5. Obtener nuevo estado (menor energía)
        perm = sortperm(Real.(evals))
        c_new = zeros(Float64, basis.num_splines)
        c_new[active] = evecs[:, perm[1]]
        
        # Normalizar nuevo estado
        c_new ./= sqrt(dot(c_new, S * c_new))

        # 6. Calcular Energía Total Correcta
        #    E_tot = 2*eps - <J> (corrección por doble conteo)
        eps = evals[perm[1]]
        E_J = dot(c_new, J_mat * c_new)
        E_total = 2 * eps - E_J

        push!(history_E, E_total)
        diff = abs(E_total - E_old)
        
        @printf("   Iter %2d: E = %.8f Ha (Δ = %.1e)\n", iter, E_total, diff)
        
        if diff < TOL
            println("   >> Convergencia alcanzada.")
            break
        end
        
        # 7. Mixing Lineal (Amortiguamiento)
        #    c_next = mix * c_old + (1-mix) * c_new
        c_current = MIXING * c_current + (1.0 - MIXING) * c_new
        
        # Re-normalizar después del mixing (CRUCIAL)
        c_current ./= sqrt(dot(c_current, S * c_current))
        
        E_old = E_total
    end

    # --------------------------------------------------------------------------
    # D. ANÁLISIS VISUAL
    # --------------------------------------------------------------------------
    println("\nD. Generando Gráficas...")
    
    r_plot = collect(range(0, R_MAX, length=300))
    
    # Usamos la nueva función eficiente evaluate_orbital
    psi_init = evaluate_orbital(basis, coeffs_initial, r_plot)
    psi_final = evaluate_orbital(basis, c_current, r_plot)
    
    # Densidades |P(r)|^2
    rho_initial = psi_init .^ 2
    rho_final   = psi_final .^ 2

    # Gráfica 1: Convergencia de Energía
    p1 = plot(history_E, marker=:circle, label="Energía SCF",
              xlabel="Iteración", ylabel="E (Ha)", title="Convergencia SCF",
              legend=:right)
    hline!(p1, [-2.86168], label="Límite HF", linestyle=:dash, color=:red)

    # Gráfica 2: Efecto de Apantallamiento
    # El orbital se "hincha" (se mueve a la derecha) debido a la repulsión
    p2 = plot(r_plot, rho_initial, label="He+ (Z_eff ≈ 2)", color=:blue, ls=:dash)
    plot!(p2, r_plot, rho_final, label="He (SCF, Z_eff < 2)", color=:red, lw=2, fill=(0, 0.1, :red))
    
    plot!(p2, title="Apantallamiento Electrónico (Expansión Radial)", 
          xlabel="Radio r (u.a.)", ylabel="Densidad Radial |P(r)|²")

    l = @layout [a; b]
    p_final = plot(p1, p2, layout=l, size=(700, 800))
    display(p_final)
    # savefig(p_final, "helio_apantallamiento_v2.png")
end

solve_helium_screening()
