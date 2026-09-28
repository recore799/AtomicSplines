# =============================================================================
#  Diagnostico: estados por debajo de -Z^2 en el Fock de core congelado del CI
#
#  `build_orbital_pool` (np2_ci_full.jl) diagonaliza, por cada l <= 3, el Fock de core congelado
#  F_l = T + V + l(l+1)/2r^2 + J_core - K_core sobre TODOS los splines de la base, y
#  `extract_virtuals` (src/ci.jl) descarta los autovalores menores que -Z^2 antes de tomar los
#  virtuales. El apendice D de la tesis habla de ese filtro; aqui se mide si alguna vez actua.
#
#  Para cada .jld2 registrado reconstruye ese mismo operador con los orbitales guardados y
#  reporta, por canal: cuantos autovalores caen por debajo de -Z^2, el menor autovalor, y el
#  mismo menor autovalor cuando se quitan el primer y el ultimo spline (los unicos no nulos en
#  r = 0 y en R_max), que es el espacio donde P(0) = P(R_max) = 0 se cumple por construccion.
#  Control: los primeros autovalores del canal deben reproducir las energias de los orbitales
#  de core de ese l, que son autovectores de su propio Fock de canal (no identico a este: el de
#  la tesis incluye el promedio de la capa abierta; la diferencia es chica y solo es un control).
#
#  No corre ningun SCF ni CI y no escribe nada. Segundos.
#  Uso: julia --project=. benchmarks/np2_sequence/tesis/diagnostico_estados_espurios.jl [--vpol]
# =============================================================================

isdefined(Main, :CABECERA_CI) || include(joinpath(@__DIR__, "comun.jl"))
isdefined(Main, :run_full_ci) || include(joinpath(NP2, "np2_ci_full.jl"))   # s_orthonormalize!

function espectro_canal(ws, core, l, idx)
    J = (build_total_J_matrix!(ws, core); copy(ws.J))
    K = zeros(size(ws.S)); assemble_K_matrix!(ws, K, l, core)
    F = ws.T .+ ws.V .+ (l * (l + 1) / 2.0) .* ws.R_inv2 .+ J .- K
    return eigen(Symmetric(F[idx, idx]), Symmetric(ws.S[idx, idx])).values
end

function main_espurios()
    man = manifiesto()
    con_vpol = "--vpol" in ARGS
    println("\n", "="^112)
    println(" Autovalores del Fock de core congelado del CI (build_orbital_pool) frente al filtro -Z^2")
    println("="^112)
    @printf("%-44s %2s %6s | %9s %14s | %14s | %s\n", "archivo", "l", "n<-Z2", "-Z^2", "min (todos)",
            "min (sin 1,N)", "core de este l: eps guardado / autovalor")
    println("-"^112)
    for (archivo, r) in sort(collect(man); by = p -> (orden_elemento(p[2]["elemento"]),
                                                      p[2]["alpha_d"], p[2]["estado"]))
        (r["estado"] == ESTADO && (r["alpha_d"] == 0.0 || con_vpol)) || continue
        el = r["elemento"]
        verificar_malla(el, archivo, man)
        ws = workspace(el; alpha_d = r["alpha_d"])
        orbs = load(ruta_np2(archivo))["orbitals"]
        core = orbs[1:end-1]
        Z = INFO[el].Z
        n = ws.basis.num_splines
        for l in 0:3
            ev_todos = espectro_canal(ws, core, l, 1:n)
            ev_int = espectro_canal(ws, core, l, 2:n-1)
            nbajo = count(<(-Z^2), ev_todos)
            cl = [o for o in core if o.l == l]
            ctrl = isempty(cl) ? "--" :
                   join([@sprintf("%.4f/%.4f", o.energy, ev_todos[nbajo + i]) for (i, o) in enumerate(cl)], " ")
            @printf("%-44s %2d %6d | %9.1f %14.4f | %14.4f | %s\n", archivo, l, nbajo, -Z^2,
                    minimum(ev_todos), minimum(ev_int), ctrl)
        end
    end
    println("-"^112)
    println(" n<-Z2: autovalores que extract_virtuals descarta.  'sin 1,N': mismo operador sin el primer")
    println(" y el ultimo spline. El control compara la energia guardada de cada orbital de core con el")
    println(" autovalor que le corresponde despues de los descartados.")
end

"""
Solapamiento de los orbitales del espacio activo con el core congelado del mismo l, con el
mismo procedimiento que `build_orbital_pool`: autovectores del Fock de core congelado sobre
todos los splines, filtro -Z^2 y salto de los `n_core(l)` primeros, sustitucion del primer
virtual de valencia por el orbital SCF y reortonormalizacion dentro del bloque. El CI supone
que esos orbitales son ortogonales al core; aqui se mide cuanto se aparta de ello.
"""
function solapamiento_con_core(ws, orbs, Z, m)
    core, val = orbs[1:end-1], orbs[end]
    n = ws.basis.num_splines
    res = Dict{Int,Tuple{Float64,Float64}}()
    for l in 0:LMAX_CI
        cl = [o for o in core if o.l == l]
        isempty(cl) && continue
        J = (build_total_J_matrix!(ws, core); copy(ws.J))
        K = zeros(size(ws.S)); assemble_K_matrix!(ws, K, l, core)
        F = ws.T .+ ws.V .+ (l * (l + 1) / 2.0) .* ws.R_inv2 .+ J .- K
        ev, vec = eigen(Symmetric(F), Symmetric(ws.S))
        block = extract_virtuals(ev, vec, l, length(cl), m, 1:n, n, ws, Z)
        l == val.l && !isempty(block) && (block[1] = deepcopy(val))
        s_orthonormalize!(block, ws)
        s = [abs(dot(v.coeffs, ws.S * c.coeffs)) for v in block, c in cl]
        # Peso total del core en el primer virtual y el maximo solapamiento individual.
        peso = maximum(sum(s .^ 2; dims = 2))
        res[l] = (maximum(s), peso)
    end
    return res
end

function main_solapamientos()
    man = manifiesto()
    println("\n", "="^86)
    println(" Solapamiento del espacio activo del CI con el core congelado del mismo l (m de produccion)")
    println("="^86)
    @printf("%-40s %4s | %s\n", "archivo", "m", "l: max |<v|c>|, max peso del core en un virtual")
    println("-"^86)
    for (archivo, r) in sort(collect(man); by = p -> (orden_elemento(p[2]["elemento"]), p[2]["alpha_d"]))
        (r["estado"] == ESTADO && r["alpha_d"] == 0.0) || continue
        el = r["elemento"]
        ws = workspace(el)
        orbs = load(ruta_np2(archivo))["orbitals"]
        m = m_produccion(el)
        s = solapamiento_con_core(ws, orbs, INFO[el].Z, m)
        @printf("%-40s %4d | %s\n", archivo, m,
                join([@sprintf("l=%d: %.2e, %.2e", l, s[l]...) for l in sort(collect(keys(s)))], "   "))
    end
    println("-"^86)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_espurios()
    main_solapamientos()
end
