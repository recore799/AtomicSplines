# =============================================================================
#  zeta efectiva que "pide" el experimento
#
#  Breit-Pauli de p^2 en acoplamiento intermedio SIN CI, con las matrices de
#  correlacion-estructura-fina.tex (sec:matrices_bp): energias LS libres y una sola zeta.
#  Tres niveles del NIST fijan las dos separaciones LS y zeta; el cuarto queda como
#  prediccion, y su error mide lo que una zeta de un cuerpo no puede describir.
#  Lo usan la etapa 3 (tabla zeta.tex) y la etapa 4 (figuras de zeta).
# =============================================================================

"""
    niveles_bp(dD, dS, zeta, D = 0)

Niveles de p^2 (cm^-1 desde 3P_0, orden [3P_0, 3P_1, 3P_2, 1D_2, 1S_0]). `D` es un termino
tensorial de rango (2,2) dentro del 3P, con el patron 5 : -5/2 : 1/2 que da
3(L.S)^2 + (3/2)(L.S) - L^2 S^2 para J = 0, 1, 2 (traza nula): la forma del espin-espin.
"""
function niveles_bp(dD, dS, zeta, D = 0.0)
    e0 = eigvals(Symmetric([(-zeta + 5D) (-sqrt(2) * zeta); (-sqrt(2) * zeta) dS]))
    e2 = eigvals(Symmetric([(zeta / 2 + D / 2) (zeta / sqrt(2)); (zeta / sqrt(2)) dD]))
    L = [e0[1], -zeta / 2 - 5D / 2, e2[1], e2[2], e0[2]]
    return L .- L[1]
end

function newton(f, x0; tol = 1e-9, max_iter = 200)
    x = float.(copy(x0))
    for _ in 1:max_iter
        r = f(x)
        norm(r) < tol && return x
        J = zeros(length(r), length(x))
        for j in eachindex(x)
            h = 1e-6 * max(1.0, abs(x[j]))
            xp = copy(x)
            xp[j] += h
            J[:, j] = (f(xp) - r) / h
        end
        x -= J \ r
    end
    error("newton: no convergio")
end

"""
    zeta_nist(el; excluir) -> (zeta, prediccion)

zeta (cm^-1) que, con energias LS libres, reproduce exactamente tres niveles del NIST.
`excluir = 2` deja fuera el 3P_1 y ajusta 3P_2, 1D_2 y 1S_0; `excluir = 3` deja fuera el
3P_2. `prediccion` es el nivel excluido segun el ajuste.
"""
function zeta_nist(el; excluir::Int = 2)
    n = NIST_NIVELES[el]
    usar = excluir == 2 ? [3, 4, 5] : [2, 4, 5]
    x = newton(x -> niveles_bp(x...)[usar] .- n[usar],
               [n[4], n[5], excluir == 2 ? 2n[3] / 3 : 2n[2]])
    return (zeta = x[3], prediccion = niveles_bp(x...)[excluir])
end

"""
    zeta_nist_tensorial(el) -> (zeta, D)

Ajuste exacto de los cuatro niveles excitados del NIST con zeta y un termino tensorial D
dentro del 3P. En carbono D mide el espin-espin; en germanio y estanio absorbe ademas lo que
una zeta de un cuerpo no describe, y no conviene leerlo como espin-espin.
"""
function zeta_nist_tensorial(el)
    n = NIST_NIVELES[el]
    x = newton(x -> niveles_bp(x...)[2:5] .- n[2:5], [n[4], n[5], 2n[3] / 3, 0.0])
    return (zeta = x[3], D = x[4])
end
