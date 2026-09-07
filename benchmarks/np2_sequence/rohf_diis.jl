# =============================================================================
#  C-DIIS por bloques para el ROHF dependiente del termino de la secuencia np^2
#
#  Adaptacion del prototipo de capa cerrada de benchmarks/closed_shell/11_radon.jl.
#  La diferencia de fondo es que alli hay UN operador de Fock por cada l, mientras
#  que aqui hay uno por cada CANAL: la capa abierta np ve su propio J/K intra-shell y
#  por eso convive con un F_core_p distinto sobre el mismo espacio activo l=1. Cada
#  canal aporta su matriz de error con su pareja (F, D) propia; mezclarlas invalida
#  el conmutador.
#
#  POR QUE EL CONMUTADOR VA PROYECTADO. El vector de error de Pulay crudo, FDS-SDF,
#  NO se anula en el punto fijo de este SCF. El orbital de valencia sale de su propio
#  Fock y despues se reortogonaliza por Gram-Schmidt contra las capas cerradas del
#  mismo l, que vienen de otro operador; en el punto fijo cumple
#
#      F_np |np> = eps |np> + sum_core lambda_c |c>,
#
#  donde los lambda_c son los multiplicadores de Lagrange de la restriccion de
#  ortogonalidad. El termino de multiplicadores deja al conmutador crudo un piso no
#  nulo (medido: ~1.2e-4 en Si), y DIIS extrapolando sobre vectores de error que ya no
#  bajan hace un paseo aleatorio en la energia en vez de converger.
#
#  Lo que si tiene que anularse es la componente VIRTUAL del gradiente orbital. Al
#  proyectar fuera del espacio ocupado del canal,
#
#      e = (1 - S C_occ C_occ^T) (F D S - S D F) = (1 - S C_occ C_occ^T) F D S,
#
#  el segundo termino muere identicamente y queda el bloque ocupado-virtual, que es
#  cero exactamente en la solucion ROHF. El conmutador crudo se sigue reportando como
#  diagnostico, porque la diferencia entre los dos es justo el argumento de la tesis.
#
#  ADVERTENCIA FISICA (va a la tesis): aun proyectado, esto es "C-DIIS por bloques
#  sobre operadores de Fock dependientes del termino", no el C-DIIS canonico de
#  Roothaan: no hay un unico Fock que diagonalizar y el acoplamiento entre la capa
#  abierta y las cerradas no entra en el residual. La convergencia se declara con la
#  ENERGIA; el residual solo decide cuando cambiar de regimen.
# =============================================================================

using LinearAlgebra

# Por encima de este condicionamiento el bloque de Pulay ya no da informacion nueva:
# los vectores de error casi se repiten y el solve amplifica el ruido.
const DIIS_COND_MAX = 1.0e12

# Coeficientes de Pulay mas grandes que esto significan cancelacion catastrofica entre
# vectores de error casi paralelos: el solve "funciona" pero la extrapolacion es ruido.
const DIIS_CMAX = 25.0

# Iteraciones seguidas sin un residual mejor antes de dar el error por estancado.
const DIIS_MAX_STALL = 6

"""
    DIISState(; max_hist, thresh)

Historial de Pulay y estado del cambio de regimen.

`thresh` es el umbral sobre max_canales ||FDS-SDF||_inf por debajo del cual se apaga
el level shift y se enciende la extrapolacion. Es un criterio por RESIDUAL, no por
numero fijo de iteraciones: el level shift es una muleta numerica y hay que soltarla
cuando el sistema ya esta en la cuenca de atraccion, no en una iteracion arbitraria.
`switch_iter` guarda en cual ocurrio, que es un dato que va a la tabla.
"""
mutable struct DIISState
    max_hist::Int
    thresh::Float64
    F_hist::Vector{Vector{Matrix{Float64}}}   # Focks CRUDAS (sin desplazar) por canal
    E_hist::Vector{Vector{Matrix{Float64}}}   # matrices de error por canal
    on::Bool
    switch_iter::Int
    resid::Float64        # residual proyectado: el que decide el cambio de regimen
    resid_raw::Float64    # conmutador crudo FDS-SDF, solo diagnostico
    restarts::Int
    frozen::Bool          # el residual toco su piso: se devuelve el control al SCF simple
    freeze_iter::Int
    best_resid::Float64
    stall::Int
end

DIISState(; max_hist::Int = 6, thresh::Float64 = 1.0e-2) =
    DIISState(max_hist, thresh, Vector{Matrix{Float64}}[], Vector{Matrix{Float64}}[],
              false, 0, Inf, Inf, 0, false, 0, Inf, 0)

"""
    block_density(orbitals, idxs, act)

Matriz de densidad del canal, restringida al espacio activo. `idxs` son los orbitales
que ocupan ese canal y entran con su ocupacion: la capa abierta np entra con la suya
fraccionaria (2 de 6), que es justo lo que la distingue del canal de core del mismo l.
"""
function block_density(orbitals, idxs, act::UnitRange{Int})
    D = zeros(Float64, length(act), length(act))
    for i in idxs
        c = orbitals[i].coeffs[act]
        D .+= orbitals[i].occ .* (c * c')
    end
    return D
end

"""
    occ_block(orbitals, idxs, act)

Matriz con TODOS los orbitales ocupados de un mismo l, restringidos al espacio activo
y puestos en columnas. Es la base del espacio fuera del cual hay que proyectar: incluye
los orbitales de los otros canales del mismo l, porque son ellos los que generan los
multiplicadores de Lagrange. Los scripts los dejan S-ortonormales por Gram-Schmidt.
"""
function occ_block(orbitals, idxs, act::UnitRange{Int})
    C = zeros(Float64, length(act), length(idxs))
    for (col, i) in enumerate(idxs)
        C[:, col] = orbitals[i].coeffs[act]
    end
    return C
end

"""
    commutator_error(F, D, S, act)

Vector de error de Pulay e = FDS - SDF sobre el espacio activo del canal. Se anula
exactamente cuando D conmuta con F en la metrica S, es decir cuando los orbitales
ocupados del canal generan un subespacio invariante de su propio Fock.
"""
function commutator_error(F::Matrix{Float64}, D::Matrix{Float64},
                          S::Matrix{Float64}, act::UnitRange{Int})
    Fa = @view F[act, act]
    Sa = @view S[act, act]
    return Fa * D * Sa - Sa * D * Fa
end

"""
    projected_error(F, D, C_occ, S, act)

Conmutador proyectado fuera del espacio ocupado del canal,
`(1 - S C_occ C_occ^T)(FDS - SDF)`. El termino `S D F` muere bajo el proyector, asi
que lo que queda es el bloque ocupado-virtual del gradiente orbital: cero exacto en la
solucion ROHF, a diferencia del conmutador crudo.
"""
function projected_error(F::Matrix{Float64}, D::Matrix{Float64}, C_occ::Matrix{Float64},
                         S::Matrix{Float64}, act::UnitRange{Int})
    Fa = @view F[act, act]
    Sa = @view S[act, act]
    M = Fa * D * Sa
    return M .- Sa * (C_occ * (C_occ' * M))
end

"""
    diis_extrapolate(st) -> Vector{Matrix} o nothing

Resuelve el sistema aumentado de Pulay

    [ B   -1 ] [ c ]   [ 0 ]
    [ -1   0 ] [ l ] = [ -1 ],   B_ij = sum_canales <e_i, e_j>_F

y devuelve F_eff = sum_i c_i F_i por canal. Devuelve `nothing` si el historial no da
para extrapolar de forma fiable; el que llama debe entonces reiniciarlo y dar una
iteracion de SCF simple. Sin esa salvaguarda DIIS diverge en silencio y sale peor que
no tenerlo.
"""
function diis_extrapolate(st::DIISState)
    while length(st.E_hist) >= 2
        m = length(st.E_hist)
        B = zeros(Float64, m + 1, m + 1)
        for i in 1:m, j in i:m
            v = 0.0
            for ch in eachindex(st.E_hist[i])
                v += dot(st.E_hist[i][ch], st.E_hist[j][ch])
            end
            B[i, j] = v
            B[j, i] = v
        end
        B[1:m, end] .= -1.0
        B[end, 1:m] .= -1.0

        # Si el bloque de Pulay se degrada, tirar los vectores mas viejos y reintentar.
        if cond(B[1:m, 1:m]) > DIIS_COND_MAX
            popfirst!(st.F_hist)
            popfirst!(st.E_hist)
            continue
        end

        rhs = zeros(Float64, m + 1)
        rhs[end] = -1.0
        c = try
            B \ rhs
        catch
            return nothing
        end

        # sum_i c_i = 1 lo impone el multiplicador de Lagrange: si no se cumple, el
        # solve no es de fiar aunque no haya lanzado excepcion.
        (all(isfinite, c) && abs(sum(@view c[1:m]) - 1.0) <= 1.0e-6) || return nothing

        # Coeficientes enormes con signos alternados = los vectores de error se volvieron
        # casi paralelos. cond(B) no siempre lo detecta, porque el ruido relativo puede
        # ser moderado aunque la extrapolacion ya no informe nada.
        if maximum(abs, @view c[1:m]) > DIIS_CMAX
            popfirst!(st.F_hist)
            popfirst!(st.E_hist)
            continue
        end

        F_eff = [zeros(Float64, size(F)) for F in st.F_hist[1]]
        for i in 1:m, ch in eachindex(F_eff)
            F_eff[ch] .+= c[i] .* st.F_hist[i][ch]
        end
        return F_eff
    end
    return nothing
end

"""
    diis_update!(st, iter, F_raw, D, C_occ, act, S) -> (F_eff, diis_on)

Un paso de C-DIIS. `F_raw` son las Focks SIN desplazar y `D[c]`, `C_occ[c]`, `act[c]`
la densidad del canal, los ocupados de su l y el espacio activo, en el mismo orden.

Mientras el residual siga por encima de `st.thresh` devuelve las Focks crudas con
`diis_on = false` y no acumula historial: el regimen de arranque (level shift + SCF
simple) queda del lado del que llama. Cuando el residual baja del umbral se enciende
la extrapolacion, y a partir de ahi el level shift ya NO debe aplicarse: sus
proyectores cambian en cada iteracion, asi que extrapolar Focks desplazadas de
iteraciones distintas seria inconsistente.
"""
function diis_update!(st::DIISState, iter::Int,
                      F_raw::Vector{Matrix{Float64}},
                      D::Vector{Matrix{Float64}},
                      C_occ::Vector{Matrix{Float64}},
                      act::Vector{UnitRange{Int}},
                      S::Matrix{Float64})

    errs = [projected_error(F_raw[c], D[c], C_occ[c], S, act[c]) for c in eachindex(F_raw)]
    st.resid = maximum(maximum(abs, e) for e in errs)
    st.resid_raw = maximum(maximum(abs, commutator_error(F_raw[c], D[c], S, act[c]))
                           for c in eachindex(F_raw))

    if !st.on
        st.resid < st.thresh || return F_raw, false
        st.on = true
        st.switch_iter = iter
    end

    # Piso del residual. El orbital de valencia se reortogonaliza por Gram-Schmidt
    # contra las capas cerradas del mismo l, y esa correccion deja al gradiente un
    # piso no nulo (medido: ~2e-7 en Si). Por debajo de el los vectores de error ya
    # no traen informacion: la extrapolacion se vuelve ruido y la energia entra en un
    # ciclo limite en vez de converger. Cuando el residual deja de mejorar se congela
    # DIIS y el tramo final lo hace el SCF simple, que si tiene el punto fijo correcto.
    if st.frozen
        return F_raw, true
    end
    if st.resid < 0.9 * st.best_resid
        st.best_resid = st.resid
        st.stall = 0
    else
        st.stall += 1
        if st.stall >= DIIS_MAX_STALL
            st.frozen = true
            st.freeze_iter = iter
            empty!(st.F_hist)
            empty!(st.E_hist)
            return F_raw, true
        end
    end

    push!(st.F_hist, [copy(F) for F in F_raw])
    push!(st.E_hist, errs)
    while length(st.F_hist) > st.max_hist
        popfirst!(st.F_hist)
        popfirst!(st.E_hist)
    end

    # Con un solo vector de error todavia no hay nada que extrapolar. Eso NO es un
    # fallo: si se tratara como tal se vaciaria el historial y nunca llegaria a dos.
    length(st.F_hist) >= 2 || return F_raw, true

    F_eff = diis_extrapolate(st)
    if F_eff === nothing
        empty!(st.F_hist)
        empty!(st.E_hist)
        st.restarts += 1
        return F_raw, true          # SCF simple, sin shift, y a reconstruir historial
    end
    return F_eff, true
end
