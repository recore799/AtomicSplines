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

# Espacio activo del CI: m orbitales por cada l = 0..LMAX_CI. El ultimo tamano de la
# curva de convergencia es el de produccion, el que va a las tablas.
const LMAX_CI = 3
const TAMANOS_CONVERGENCIA = [4, 8, 12, 16, 20]
const M_PRODUCCION = last(TAMANOS_CONVERGENCIA)
espacio_activo(m::Int) = Dict(l => m for l in 0:LMAX_CI)

# --- Literatura -----------------------------------------------------------------

# NIST ASD: niveles de C I, Si I, Ge I y Sn I en cm^-1 referidos a 3P_0, en el orden
# [3P_0, 3P_1, 3P_2, 1D_2, 1S_0]. Son los de generate_table.jl, redondeados a 0.1 cm^-1.
# PENDIENTE: la tesis tiene que citar la version del ASD y la fecha de consulta.
const NIST_NIVELES = Dict(
    "C"  => [0.0,   16.4,   43.4, 10192.6, 21648.0],
    "Si" => [0.0,   77.1,  223.2,  6298.8, 15394.4],
    "Ge" => [0.0,  557.1, 1410.0,  7125.3, 16367.1],
    "Sn" => [0.0, 1691.8, 3427.7,  8613.0, 17162.6],
)

# Primer potencial de ionizacion (eV). C, Si y Ge son los de tab:koopmans; el del Sn es el
# valor estandar y hay que confirmarlo en el NIST ASD antes de publicarlo.
const NIST_IONIZACION = Dict("C" => 11.26, "Si" => 8.15, "Ge" => 7.90, "Sn" => 7.34)

# Factor g de Lande del 3P_2. Solo el de germanio aparece en la tesis (discusion.tex:67),
# y hay que VERIFICARLO; el NIST ASD da los cuatro.
const NIST_LANDE_3P2 = Dict{String,Union{Float64,Nothing}}(
    "C" => nothing, "Si" => nothing, "Ge" => 1.496, "Sn" => nothing)

# Froese Fischer (1977), limite Hartree-Fock del termino 3P, copiado tal cual de
# resultados.tex (tab:resultados_energia_global, tab:momentos_inversos y
# tab:integrales_slater). Van como texto para conservar las cifras del libro; `nothing`
# donde no hay un dato confiable:
#   - F2(4p,4p) del Ge: la tesis imprime 0.16593265, que es copia del de Si. Recuperarlo
#     del libro.
#   - Virial del Sn: la tesis pone 2.00000002 en la columna Delta, quiza la referencia
#     desplazada. Confirmar en el libro.
const FF = Dict(
    "C"  => (E = "-37.688619", T = "37.688619", virial = "1.999999998",
             r3 = "1.69181", F0 = "0.53860360", F2 = "0.24330170"),
    "Si" => (E = "-288.85436", T = "288.85437", virial = "1.999999965",
             r3 = "2.05423", F0 = "0.32971786", F2 = "0.16593265"),
    "Ge" => (E = "-2075.3597", T = "2075.3597", virial = "2.00000005",
             r3 = "4.80120", F0 = "0.32624597", F2 = nothing),
    "Sn" => (E = "-6022.9317", T = "6022.9316", virial = nothing,
             r3 = "6.83248", F0 = "0.27875374", F2 = "0.14722241"),
)
