# Valores citados en el texto

GENERADO por `benchmarks/np2_sequence/tesis/etapa3_tablas.jl`; no editar a mano.
Convencion de orbitales: 3P (sensibilidad: av); m = 20, lmax = 3.

## Pozo del potencial radial efectivo V_rad = V_eff + 1/r^2

- Carbono: minimo -1.02 Ha en r = 0.490 a0
- Silicio: minimo -15.77 Ha en r = 0.165 a0
- Germanio: minimo -138.88 Ha en r = 0.067 a0
- Estaño: minimo -399.43 Ha en r = 0.042 a0

## zeta y F^2 (cm^-1)

- Carbono: zeta(3P) = 42.216, zeta(av) = 41.480 (razon 1.0178), F^2 = 53398.6
- Silicio: zeta(3P) = 140.992, zeta(av) = 138.903 (razon 1.0150), F^2 = 36421.4
- Germanio: zeta(3P) = 823.671, zeta(av) = 810.736 (razon 1.0160), F^2 = 35713.9
- Estaño: zeta(3P) = 1885.936, zeta(av) = 1858.739 (razon 1.0146), F^2 = 32313.7

## zeta que pide el NIST con una sola zeta (Breit-Pauli p^2 sin CI)

- Carbono: sin 3P_1 -> zeta = 28.9 (predice 3P_1 = 14.5, NIST 16.4); sin 3P_2 -> zeta = 32.6 (predice 3P_2 = 49.0, NIST 43.4)
- Silicio: sin 3P_1 -> zeta = 148.1 (predice 3P_1 = 76.9, NIST 77.1); sin 3P_2 -> zeta = 148.5 (predice 3P_2 = 223.8, NIST 223.2)
- Germanio: sin 3P_1 -> zeta = 920.6 (predice 3P_1 = 564.5, NIST 557.1); sin 3P_2 -> zeta = 910.4 (predice 3P_2 = 1394.2, NIST 1410.0)
- Estaño: sin 3P_1 -> zeta = 2241.1 (predice 3P_1 = 1727.3, NIST 1691.8); sin 3P_2 -> zeta = 2207.5 (predice 3P_2 = 3382.4, NIST 3427.7)

## zeta con termino tensorial dentro del 3P (ajuste exacto de los cuatro niveles)

- Carbono: zeta = 28.0, D = -0.31
- Silicio: zeta = 148.0, D = -0.03
- Germanio: zeta = 924.6, D = 1.37
- Estaño: zeta = 2265.4, D = 8.19

## CI de valencia (alpha_d = 0)

- Carbono 3P: E_corr = -8.698045e-03 Ha (-1909.0 cm^-1), g_eff(3P_2) = 1.499996, niveles = 0.00 / 21.15 / 63.08 / 11148.39 / 24574.78
- Carbono av: E_corr = -8.965971e-03 Ha (-1967.8 cm^-1), g_eff(3P_2) = 1.499997, niveles = 0.00 / 20.77 / 61.96 / 11127.33 / 24524.70
- Silicio 3P: E_corr = -8.562860e-03 Ha (-1879.3 cm^-1), g_eff(3P_2) = 1.499903, niveles = 0.00 / 72.35 / 210.58 / 7285.25 / 14674.51
- Silicio av: E_corr = -8.722626e-03 Ha (-1914.4 cm^-1), g_eff(3P_2) = 1.499906, niveles = 0.00 / 71.21 / 207.32 / 7276.55 / 14658.14
- Germanio 3P: E_corr = -7.060409e-03 Ha (-1549.6 cm^-1), g_eff(3P_2) = 1.496164, niveles = 0.00 / 490.74 / 1257.58 / 7851.64 / 15780.18
- Germanio av: E_corr = -7.219809e-03 Ha (-1584.6 cm^-1), g_eff(3P_2) = 1.496278, niveles = 0.00 / 481.99 / 1237.04 / 7828.27 / 15748.30
- Estaño 3P: E_corr = -6.244565e-03 Ha (-1370.5 cm^-1), g_eff(3P_2) = 1.472668, niveles = 0.00 / 1367.10 / 2921.70 / 8710.14 / 16045.30
- Estaño av: E_corr = -6.379331e-03 Ha (-1400.1 cm^-1), g_eff(3P_2) = 1.473394, niveles = 0.00 / 1343.57 / 2878.98 / 8659.29 / 15987.93

## Barridos de V_pol: error frente al NIST

- Germanio alpha_d = 0.25: 3P_1 -7.7%, 3P_2 -6.9%, 1D_2 = 7967.9, 1S_0 = 15999.0
- Germanio alpha_d = 0.50: 3P_1 -3.2%, 3P_2 -2.8%, 1D_2 = 8089.1, 1S_0 = 16226.7
- Germanio alpha_d = 0.75: 3P_1 +1.6%, 3P_2 +1.5%, 1D_2 = 8215.3, 1S_0 = 16463.5  <- elegido
- Germanio alpha_d = 1.00: 3P_1 +6.6%, 3P_2 +6.0%, 1D_2 = 8346.7, 1S_0 = 16709.7
- Estaño alpha_d = 1.00: 3P_1 -14.6%, 3P_2 -10.5%, 1D_2 = 8980.8, 1S_0 = 16457.5
- Estaño alpha_d = 2.00: 3P_1 -8.8%, 3P_2 -5.4%, 1D_2 = 9301.8, 1S_0 = 16952.1
- Estaño alpha_d = 3.00: 3P_1 -1.9%, 3P_2 +0.8%, 1D_2 = 9675.4, 1S_0 = 17530.5  <- elegido
- Estaño alpha_d = 4.00: 3P_1 +6.4%, 3P_2 +8.0%, 1D_2 = 10104.5, 1S_0 = 18195.5
- Estaño alpha_d = 5.00: 3P_1 +16.2%, 3P_2 +16.3%, 1D_2 = 10592.4, 1S_0 = 18950.1
- Estaño alpha_d = 6.00: 3P_1 +27.5%, 3P_2 +25.8%, 1D_2 = 11142.0, 1S_0 = 19797.7
