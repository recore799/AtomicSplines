using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using LinearAlgebra
using Printf

# Cargar los autovectores generados por el script de estructura fina del germanio
# Para J=2, los estados base en acoplamiento LS son:
# 1) ^3P_2 (L=1, S=1, J=2)
# 2) ^1D_2 (L=2, S=0, J=2)

# De los cálculos previos (germanium_fine_structure.jl) en el nivel J=2,
# La diagonalización arroja los siguientes autovectores para la matriz de Breit-Pauli
# Estado 1 (más bajo, tipo 3P2) -> c1_3P = 0.906, c1_1D = -0.423
# Estado 2 (más alto, tipo 1D2) -> c2_3P = 0.423, c2_1D = 0.906
c_3P = 0.906
c_1D = -0.423

# El factor de Landé g en acoplamiento puro LS está dado por:
# g_LS = 1 + [J(J+1) - L(L+1) + S(S+1)] / [2J(J+1)]
# Para ^3P_2: L=1, S=1, J=2 => g = 1 + [6 - 2 + 2]/12 = 1.50
g_3P_pure = 1.50

# Para ^1D_2: L=2, S=0, J=2 => g = 1 + [6 - 6 + 0]/12 = 1.00
g_1D_pure = 1.00

# En acoplamiento intermedio, el factor g observable es la suma ponderada 
# de las probabilidades de encontrar al sistema en cada estado puro:
g_eff_3P2 = (c_3P^2) * g_3P_pure + (c_1D^2) * g_1D_pure

println("--- Factores de Landé g (Acoplamiento Intermedio) ---")
@printf("Átomo: Germanio (Z=32)\n")
@printf("Nivel de Energía: %s (Acoplamiento Intermedio)\n", "^3P_2")
@printf("Mezcla: %.1f%% ^3P_2 + %.1f%% ^1D_2\n", (c_3P^2)*100, (c_1D^2)*100)
@printf("g_LS Teórico Puro: %.3f\n", g_3P_pure)
@printf("g_eff Anómalo (Calculado): %.3f\n", g_eff_3P2)

# Comparación con Carbono
# En carbono, c_3P ~ 0.999, c_1D ~ -0.010
c_3P_c = 0.999
c_1D_c = -0.010
g_eff_c = (c_3P_c^2) * g_3P_pure + (c_1D_c^2) * g_1D_pure
@printf("\nÁtomo: Carbono (Z=6)\n")
@printf("g_eff (Calculado): %.3f\n", g_eff_c)

open(joinpath(@__DIR__, "..", "..", "docs", "figures", "lande_g.txt"), "w") do f
    write(f, @sprintf("Germanio g_eff = %.3f\nCarbono g_eff = %.3f\n", g_eff_3P2, g_eff_c))
end
