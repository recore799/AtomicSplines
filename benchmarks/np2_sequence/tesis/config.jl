# =============================================================================
#  Configuracion de la regeneracion de resultados de la tesis
#
#  UNICO lugar donde viven las decisiones (que orbitales, que alpha_d, que espacio
#  activo) y los datos de literatura que el codigo NO puede calcular (NIST, Froese
#  Fischer). Todo lo demas se calcula de los .jld2 registrados en benchmarks/RESULTS.toml.
#
#  Despues de editar este archivo, correr la etapa 1 sin --correr: lista que SCF hacen
#  falta para la configuracion nueva y cuanto tardan, sin calcular nada.
# =============================================================================

# --- Decisiones -----------------------------------------------------------------

# Orbitales sobre los que se construye todo: CI, Breit-Pauli, V_pol, ionizacion y
# figuras. "3P" son los optimizados para el termino fundamental, los mismos que se
# contrastan con Froese Fischer, de modo que la tesis entera usa una sola convencion.
# El promedio de configuracion ("av") entra solo como prueba de sensibilidad a
# alpha_d = 0. Ver docs/claude/PLAN-REGENERACION.md, seccion 2.1.
const ESTADO = "3P"
const ESTADO_SENSIBILIDAD = "av"

const R_MAX = 30.0

# Energia Hartree-Fock que reportan las tablas. :autovalor es la que guardan los .jld2,
# E = sum_i (occ_i/2)(h_ii + eps_i). :rayleigh sustituye eps_i por <c_i|F|c_i>, que es lo
# consistente cuando Gram-Schmidt deja al orbital de valencia fuera de su autovector. Las
# dos se evaluan sobre los MISMOS orbitales, sin correr ningun SCF: ver
# diagnostico_rayleigh.jl y docs/claude/PLAN-REGENERACION.md, seccion 2.2.
const ENERGIA_HF = :autovalor

# Barridos del potencial de polarizacion del core, con r_c = 1.0 a0. Germanio: el 0.5 de
# la tesis quedaba en el borde de su barrido anterior, y 0.25-1.0 lo encierra. Estanio: los
# puntos 2-6 de tab:estanio_ci mas el 1, que ya existia en 3P. Ver seccion 2.3 del plan.
const R_C = 1.0
const BARRIDO_ALPHA = Dict("Ge" => [0.25, 0.5, 0.75, 1.0],
                           "Sn" => [1.0, 2.0, 3.0, 4.0, 5.0, 6.0])

# Criterio para el alpha_d de la columna CI+V_pol. :rms_3P minimiza la raiz del error
# relativo cuadratico medio de 3P_1 y 3P_2 frente al NIST; :solo_3P2 usa solo el 3P_2.
const CRITERIO_ALPHA = :rms_3P

# Espacio activo del CI: m orbitales por cada l = 0..LMAX_CI. El ultimo tamano de la curva
# de cada elemento es el de produccion, el que va a las tablas. El carbono converge mucho mas
# lento que el resto: a m = 20 sus singletes todavia bajaban ~200 (1D_2) y ~530 cm^-1 (1S_0)
# por paso (docs/claude/REVISION-RESULTADOS-2026-09-13.md, seccion 2.2), por eso su curva sigue.
# Costo del carbono: las R^k crecen como m^4 (~39 millones a m = 40) y la cache usa ~17 bytes
# por casilla, del orden de 1.1 GB y hasta ~1.7 GB al redimensionarse. La curva completa tarda
# ~35 min; cada tamano se guarda al terminar, asi que se puede cortar antes.
const LMAX_CI = 3
const TAMANOS_CONVERGENCIA = Dict("C"  => [4, 8, 12, 16, 20, 24, 28, 32, 36, 40],
                                  "Si" => [4, 8, 12, 16, 20],
                                  "Ge" => [4, 8, 12, 16, 20],
                                  "Sn" => [4, 8, 12, 16, 20])
m_produccion(el) = last(TAMANOS_CONVERGENCIA[el])
# La prueba de sensibilidad (promedio de configuracion frente a ESTADO) se compara al mismo m
# en los cuatro elementos.
const M_SENSIBILIDAD = 20
espacio_activo(m::Int) = Dict(l => m for l in 0:LMAX_CI)

# --- Literatura -----------------------------------------------------------------

# NIST Atomic Spectra Database, ver. 5.12 (Kramida, Ralchenko, Reader y NIST ASD Team, 2024;
# DOI 10.18434/T4W30F), consultada el 2026-09-13 en https://physics.nist.gov/asd. Niveles de
# C I, Si I, Ge I y Sn I en cm^-1 referidos a 3P_0, en el orden [3P_0, 3P_1, 3P_2, 1D_2, 1S_0],
# con todas las cifras que da la base. Fuentes que lista el ASD: C I, Haris y Kramida (2017);
# Si I, Martin y Zalubas (1983); Ge I, Sugar y Musgrove (1993); Sn I, Brill (1964) y Brown et
# al. (1977). Antes del 2026-09-13 aqui habia valores a 0.1 cm^-1 copiados de generate_table.jl;
# el 1D_2 del carbono y el 1S_0 de germanio y estanio no eran el redondeo del valor de la base
# (hasta 0.23 cm^-1 de diferencia).
const NIST_NIVELES = Dict(
    "C"  => [0.0,   16.4167130,   43.4134567, 10192.657,  21648.030],
    "Si" => [0.0,   77.115,      223.157,      6298.850,  15394.370],
    "Ge" => [0.0,  557.1341,    1409.9609,     7125.2989, 16367.3332],
    "Sn" => [0.0, 1691.806,     3427.673,      8612.955,  17162.499],
)

# Primer potencial de ionizacion (eV), del 3P_0 del neutro al 2P_1/2 del ion, NIST ASD ver. 5.12
# (misma consulta). Incertidumbres: C 1.1e-6, Si 3e-5, Ge 1.2e-5 y Sn 1.2e-5 eV.
const NIST_IONIZACION = Dict("C" => 11.2602880, "Si" => 8.15168, "Ge" => 7.899435, "Sn" => 7.343918)

# Factor g de Lande del 3P_2, NIST ASD ver. 5.12 (misma consulta). Son medidos, con el g_s real
# del electron, y por eso la etapa 3 los compara con g calculado con g_s = G_S. Si I no tiene
# factores g en el ASD. C I: 1.5010469(50). Ge I: 1.49458, sin incertidumbre, bajo la fuente
# primaria que lista el ASD (Sugar y Musgrove, 1993). Sn I: 1.452, de Moore (1958). La tesis
# citaba 1.496 para el germanio (discusion.tex:67), que no es el valor del ASD.
const NIST_LANDE_3P2 = Dict{String,Union{Float64,Nothing}}(
    "C" => 1.5010469, "Si" => nothing, "Ge" => 1.49458, "Sn" => 1.452)

# Factor g del espin del electron libre, en valor absoluto (CODATA 2018).
const G_S = 2.00231930436256

# Froese Fischer (1977), limite Hartree-Fock del termino 3P. Van como texto para conservar las
# cifras del libro. Casi todo se copio de resultados.tex (tab:resultados_energia_global,
# tab:momentos_inversos y tab:integrales_slater); el usuario tomo del libro el 2026-09-13:
#   - F2(4p,4p) del Ge = 0.16271717 (la tesis imprimia 0.16593265, copia del de Si);
#   - virial del Sn = 2.0000000019 (la tesis lo tenia corrido a la columna Delta);
#   - zeta(2p) del C = 31.946 cm^-1 (sin confirmar si es el de Blume-Watson).
# POR VERIFICAR en el libro: el F0(4p,4p) del Ge (0.32624597) difiere 9.7e-3 del calculado,
# cuando su F2 coincide a 7e-6 y los F0 de Si y Sn a ~1e-5, asi que probablemente es otra
# transcripcion equivocada. Los <r^-3> de Si, Ge y Sn difieren 0.2-0.3 %; puede ser la malla
# cerca del nucleo o transcripcion.
const FF = Dict(
    "C"  => (E = "-37.688619", T = "37.688619", virial = "1.999999998",
             r3 = "1.69181", F0 = "0.53860360", F2 = "0.24330170", zeta = "31.946"),
    "Si" => (E = "-288.85436", T = "288.85437", virial = "1.999999965",
             r3 = "2.05423", F0 = "0.32971786", F2 = "0.16593265", zeta = nothing),
    "Ge" => (E = "-2075.3597", T = "2075.3597", virial = "2.00000005",
             r3 = "4.80120", F0 = "0.32624597", F2 = "0.16271717", zeta = nothing),
    "Sn" => (E = "-6022.9317", T = "6022.9316", virial = "2.0000000019",
             r3 = "6.83248", F0 = "0.27875374", F2 = "0.14722241", zeta = nothing),
)
