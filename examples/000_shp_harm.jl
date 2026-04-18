using Plots
gr() # Usamos el backend GR que es excelente para renderizado rápido y exportación a PDF

# Función principal para calcular y generar el objeto de superficie de un armónico
function generar_armonico(l, m)
    # Definición de la malla angular (resolución de la superficie)
    θ = range(0, π, length=100)
    φ = range(0, 2π, length=100)
    
    # Creación de matrices 2D para evaluar las funciones en toda la malla
    Θ = [t for t in θ, p in φ]
    Φ = [p for t in θ, p in φ]
    
    # Evaluación de la parte angular de los orbitales reales
    if l == 0 && m == 0
        # Orbital s
        R = ones(size(Θ))
        
    elseif l == 1
        # Orbitales p
        if m == -1
            R = sin.(Θ) .* sin.(Φ) # p_y
        elseif m == 0
            R = cos.(Θ)            # p_z
        elseif m == 1
            R = sin.(Θ) .* cos.(Φ) # p_x
        end
        
    elseif l == 2
        # Orbitales d
        if m == -2
            R = sin.(Θ).^2 .* sin.(2 .* Φ) # d_xy
        elseif m == -1
            R = sin.(2 .* Θ) .* sin.(Φ)    # d_yz
        elseif m == 0
            R = 3 .* cos.(Θ).^2 .- 1       # d_z2
        elseif m == 1
            R = sin.(2 .* Θ) .* cos.(Φ)    # d_xz
        elseif m == 2
            R = sin.(Θ).^2 .* cos.(2 .* Φ) # d_x2-y2
        end

    elseif l == 3
        # Orbitales f (l = 3, 7 orientaciones espaciales)
        if m == 0
            # Orbital f_z3
            R = 5 .* cos.(Θ).^3 .- 3 .* cos.(Θ)
        elseif m == 1
            # Orbital f_xz2
            R = sin.(Θ) .* (5 .* cos.(Θ).^2 .- 1) .* cos.(Φ)
        elseif m == -1
            # Orbital f_yz2
            R = sin.(Θ) .* (5 .* cos.(Θ).^2 .- 1) .* sin.(Φ)
        elseif m == 2
            # Orbital f_z(x2-y2)
            R = sin.(Θ).^2 .* cos.(Θ) .* cos.(2 .* Φ)
        elseif m == -2
            # Orbital f_xyz
            R = sin.(Θ).^2 .* cos.(Θ) .* sin.(2 .* Φ)
        elseif m == 3
            # Orbital f_x(x2-3y2)
            R = sin.(Θ).^3 .* cos.(3 .* Φ)
        elseif m == -3
            # Orbital f_y(3x2-y2)
            R = sin.(Θ).^3 .* sin.(3 .* Φ)
        end

    else
        R = zeros(size(Θ))
    end
    
    # La magnitud de la función define la distancia desde el origen
    r = abs.(R)
    
    # Transformación a coordenadas cartesianas
    X = r .* sin.(Θ) .* cos.(Φ)
    Y = r .* sin.(Θ) .* sin.(Φ)
    Z = r .* cos.(Θ)
    
    # El color representará el signo de la fase de la función de onda (+ o -)
    C = sign.(R)
    
    # Generación del gráfico de superficie 3D
    surface(X, Y, Z, 
            surfacecolor=C, 
            camera=(45, 30),     # Ángulo de visión óptimo
            legend=false, 
            axis=false, 
            grid=false, 
            showaxis=false,
            title="m_l = $m",
            titlefontsize=10)
end

# -------------------------------------------------------------------
# Generación de los archivos PDF solicitados
# -------------------------------------------------------------------

println("Generando el archivo para los orbitales s (l=0)...")
# l = 0 (s) - 1 orientación
p_s = generar_armonico(0, 0)
# Envolvemos el gráfico en un layout para añadir un título general
plot_final_s = plot(p_s, layout=(1, 1), size=(600, 500), plot_title="Orbital s (l = 0)")
savefig(plot_final_s, "orbitales_s.pdf")

println("Generando el archivo para los orbitales p (l=1)...")
# l = 1 (p) - 3 orientaciones
p_p_minus1 = generar_armonico(1, -1)
p_p_0      = generar_armonico(1, 0)
p_p_plus1  = generar_armonico(1, 1)

plot_final_p = plot(p_p_minus1, p_p_0, p_p_plus1, 
                    layout=(1, 3), 
                    size=(1500, 500), 
                    plot_title="Orbitales p (l = 1)")
savefig(plot_final_p, "orbitales_p.pdf")

println("Generando el archivo para los orbitales d (l=2)...")
# l = 2 (d) - 5 orientaciones
p_d_minus2 = generar_armonico(2, -2)
p_d_minus1 = generar_armonico(2, -1)
p_d_0      = generar_armonico(2, 0)
p_d_plus1  = generar_armonico(2, 1)
p_d_plus2  = generar_armonico(2, 2)

plot_final_d = plot(p_d_minus2, p_d_minus1, p_d_0, p_d_plus1, p_d_plus2, 
                    layout=(1, 5), 
                    size=(2500, 500), 
                    plot_title="Orbitales d (l = 2)")
savefig(plot_final_d, "orbitales_d.pdf")

println("Generando el archivo para los orbitales f (l=3)...")
# l = 3 (f) - 7 orientaciones
p_f_minus3 = generar_armonico(3, -3)
p_f_minus2 = generar_armonico(3, -2)
p_f_minus1 = generar_armonico(3, -1)
p_f_0      = generar_armonico(3, 0)
p_f_plus1  = generar_armonico(3, 1)
p_f_plus2  = generar_armonico(3, 2)
p_f_plus3  = generar_armonico(3, 3)

plot_final_f = plot(p_f_minus3, p_f_minus2, p_f_minus1, p_f_0, p_f_plus1, p_f_plus2, p_f_plus3,
                    layout=(1, 7), 
                    size=(3500, 500), 
                    plot_title="Orbitales f (l = 3)")
savefig(plot_final_f, "orbitales_f.pdf")






println("¡Archivos PDF generados con éxito, Sir!")
