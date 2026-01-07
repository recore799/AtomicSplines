using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AtomicSplines
using LinearAlgebra
using Plots
using Printf

function solve_helium_exploratory()
    println("\n================================================================")
    println("   EXPERIMENTO: HELIO Y EL EFECTO DE APANTALLAMIENTO")
    println("================================================================\n")

    # PARÁMETROS
    Z_ATOM = 2.0
    R_MAX = 10.0      
    N_ELEMS = 80      
    ORDER = 5
    GAMMA = 2.0

    # GENERACIÓN DE BASE
    basis = generate_basis(R_MAX, N_ELEMS, ORDER, γ=GAMMA)
    println("A. > Base: $N_ELEMS splines (Orden $ORDER) en caja de $R_MAX u.a.")

    # OPERADORES NUCLEARES
    T, V, S = assemble_core(basis, Z_ATOM)
    H_core = T + V 
    free_dofs = 2:(basis.num_splines - 1)

    # --------------------------------------------------------------------------
    # GUESS INICIAL (Átomo Hidrogenoide He+)
    # Ignoramos la repulsión e-e. Es como si los electrones se ignoracen.
    # --------------------------------------------------------------------------
    println("\nB. Calculando Estado Inicial (Sin interacción)...")
    
    vals, vecs = eigen(Matrix(H_core[free_dofs, free_dofs]), Matrix(S[free_dofs, free_dofs]))
    
    c_initial = zeros(basis.num_splines)
    c_initial[free_dofs] = vecs[:, 1]
    
    # Normalizar
    norm = sqrt(dot(c_initial, S * c_initial))
    c_initial ./= norm
    
    # Guardamos esta energía para comparar
    E_hydrogenic = 2.0 * vals[1] # 2 electrones x Energía del 1s de He+
    println("   > Energía He+ (x2): $E_hydrogenic Ha (Muy baja, ignora repulsión)")

    # --------------------------------------------------------------------------
    # CICLO SCF (Self-Consistent Field)
    # Encendemos la repulsión y dejamos que el sistema se relaje.
    # --------------------------------------------------------------------------
    println("\nC. Iniciando SCF (Encendiendo repulsión)...")
    
    c_current = copy(c_initial)
    MAX_ITER = 30
    MIXING = 0.6
    
    history_E = Float64[] # Para guardar el historial y graficar
    E_old = E_hydrogenic
    
    for iter in 1:MAX_ITER
        # Poisson: Cada electrón ve la nube de carga del otro
        y_coeffs = solve_poisson_potential(basis, c_current, T)
        J_mat = assemble_J_matrix(basis, y_coeffs)
        
        # Fock: H_eff = H_nuc + J
        F_mat = H_core + J_mat
        
        # Diagonalizar
        evals, evecs = eigen(Matrix(F_mat[free_dofs, free_dofs]), Matrix(S[free_dofs, free_dofs]))
        
        # Actualizar Coeficientes
        c_new_inner = evecs[:, 1]
        c_new = zeros(basis.num_splines)
        c_new[free_dofs] = c_new_inner
        c_new ./= sqrt(dot(c_new, S * c_new)) # Normalizar
        
        # Calcular Energía Total
        epsilon = evals[1]
        E_J = dot(c_new, J_mat * c_new) # Energía de repulsión
        E_total = 2 * epsilon - E_J
        
        push!(history_E, E_total)
        
        diff = abs(E_total - E_old)
        @printf("   Iter %2d: E = %.6f Ha (Δ = %.1e)\n", iter, E_total, diff)
        
        if diff < 1e-8
            println("   >> Convergencia alcanzada.")
            break
        end
        
        # Mezcla para estabilidad
        c_current = MIXING * c_new + (1 - MIXING) * c_current
        c_current ./= sqrt(dot(c_current, S * c_current))
        E_old = E_total
    end

    # --------------------------------------------------------------------------
    # ANÁLISIS VISUAL
    # --------------------------------------------------------------------------
    println("\nD. Generando Gráficas de Análisis...")
    
    r_plot = range(0, 5.0, length=200) # Zoom en la zona cercana al núcleo
    
    # Reconstruir densidades radiales |P(r)|²
    rho_initial = [eval_expansion(c_initial, basis, r)^2 for r in r_plot]
    rho_final   = [eval_expansion(c_current, basis, r)^2 for r in r_plot]

    # GRÁFICA 1: Convergencia de Energía
    p1 = plot(1:length(history_E), history_E, 
              marker=:circle, label="Energía SCF",
              xlabel="Iteración", ylabel="Energía (Ha)",
              title="Convergencia del Método SCF")
    # Línea de referencia Hartree-Fock
    hline!(p1, [-2.86168], label="Límite HF", linestyle=:dash, color=:red)

    # GRÁFICA 2: Efecto de Apantallamiento (La física interesante)
    p2 = plot(r_plot, rho_initial, 
              label="He+ (Sin Interacción)", 
              color=:blue, linestyle=:dash, lw=2,
              fill=(0, 0.1, :blue))
              
    plot!(p2, r_plot, rho_final, 
          label="Helio (SCF Convergido)", 
          color=:red, lw=2,
          fill=(0, 0.1, :red))
          
    plot!(p2, title="Efecto de Apantallamiento",
          xlabel="Distancia r (u.a.)", ylabel="Densidad Radial |P(r)|²")

    # Layout final
    l = @layout [a; b]
    final_plot = plot(p1, p2, layout=l, size=(700, 800))
    display(final_plot)
    
end

solve_helium_exploratory()
