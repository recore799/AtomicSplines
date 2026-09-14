# Valores citados en el texto

GENERADO por `benchmarks/np2_sequence/tesis/etapa3_tablas.jl`; no editar a mano.
Convencion de orbitales: 3P (sensibilidad: av a m = 20); lmax = 3; m de produccion: C 40, Si 20, Ge 20, Sn 20.

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
- Germanio: sin 3P_1 -> zeta = 920.6 (predice 3P_1 = 564.5, NIST 557.1); sin 3P_2 -> zeta = 910.4 (predice 3P_2 = 1394.3, NIST 1410.0)
- Estaño: sin 3P_1 -> zeta = 2241.1 (predice 3P_1 = 1727.3, NIST 1691.8); sin 3P_2 -> zeta = 2207.5 (predice 3P_2 = 3382.4, NIST 3427.7)

## zeta con termino tensorial dentro del 3P (ajuste exacto de los cuatro niveles)

- Carbono: zeta = 28.0, D = -0.31
- Silicio: zeta = 148.0, D = -0.04
- Germanio: zeta = 924.5, D = 1.37
- Estaño: zeta = 2265.4, D = 8.18

## CI de valencia (alpha_d = 0, m de produccion)

- Carbono 3P m = 40: E_corr = -9.851678e-03 Ha (-2162.2 cm^-1), g(3P_2) con g_s = 2: 1.499996, con g_s real: 1.501156, niveles = 0.00 / 21.15 / 63.08 / 10810.12 / 23894.03
- Silicio 3P m = 20: E_corr = -8.562860e-03 Ha (-1879.3 cm^-1), g(3P_2) con g_s = 2: 1.499903, con g_s real: 1.501063, niveles = 0.00 / 72.35 / 210.58 / 7285.25 / 14674.51
- Germanio 3P m = 20: E_corr = -7.060409e-03 Ha (-1549.6 cm^-1), g(3P_2) con g_s = 2: 1.496164, con g_s real: 1.497315, niveles = 0.00 / 490.74 / 1257.58 / 7851.64 / 15780.18
- Estaño 3P m = 20: E_corr = -6.244565e-03 Ha (-1370.5 cm^-1), g(3P_2) con g_s = 2: 1.472668, con g_s real: 1.473764, niveles = 0.00 / 1367.10 / 2921.70 / 8710.14 / 16045.30

## Sensibilidad a la convencion de orbitales (alpha_d = 0, m = 20)

- Carbono: av frente a 3P: 3P_1 -1.78%, 3P_2 -1.77%, 1D_2 -0.19%, 1S_0 -0.20%
- Silicio: av frente a 3P: 3P_1 -1.58%, 3P_2 -1.55%, 1D_2 -0.12%, 1S_0 -0.11%
- Germanio: av frente a 3P: 3P_1 -1.78%, 3P_2 -1.63%, 1D_2 -0.30%, 1S_0 -0.20%
- Estaño: av frente a 3P: 3P_1 -1.72%, 3P_2 -1.46%, 1D_2 -0.58%, 1S_0 -0.36%

## Convergencia de los singletes con el espacio activo (3P, alpha_d = 0)

Con los tres ultimos tamanos de cada curva: extrapolacion geometrica (incrementos en razon constante) y de potencia (E(m) = E_inf + A m^-p, ajuste exacto). Una razon cercana a 1 hace inservible la geometrica.

- Carbono 1D_2: m = 32, 36, 40 -> 10864.7, 10831.1, 10810.1 cm^-1 | geometrica 10775.4 (+5.7%, razon 0.62) | potencia 10758.2 (+5.5%, p = 3.2) | NIST 10192.7
- Carbono 1S_0: m = 32, 36, 40 -> 23988.4, 23929.9, 23894.0 cm^-1 | geometrica 23836.8 (+10.1%, razon 0.61) | potencia 23809.4 (+10.0%, p = 3.4) | NIST 21648.0
- Silicio 1D_2: m = 12, 16, 20 -> 7491.1, 7340.1, 7285.2 cm^-1 | geometrica 7254.0 (+15.2%, razon 0.36) | potencia 7225.8 (+14.7%, p = 2.9) | NIST 6298.9
- Silicio 1S_0: m = 12, 16, 20 -> 14944.7, 14767.5, 14674.5 cm^-1 | geometrica 14571.6 (-5.3%, razon 0.53) | potencia 14443.0 (-6.2%, p = 1.5) | NIST 15394.4
- Germanio 1D_2: m = 12, 16, 20 -> 8021.9, 7893.3, 7851.6 cm^-1 | geometrica 7831.7 (+9.9%, razon 0.32) | potencia 7814.3 (+9.7%, p = 3.4) | NIST 7125.3
- Germanio 1S_0: m = 12, 16, 20 -> 16069.2, 15849.5, 15780.2 cm^-1 | geometrica 15748.2 (-3.8%, razon 0.32) | potencia 15720.6 (-4.0%, p = 3.5) | NIST 16367.3
- Estaño 1D_2: m = 12, 16, 20 -> 8764.2, 8736.1, 8710.1 cm^-1 | geometrica 8388.8 (-2.6%, razon 0.93) | potencia no aplica | NIST 8613.0
- Estaño 1S_0: m = 12, 16, 20 -> 16181.5, 16090.6, 16045.3 cm^-1 | geometrica 16000.4 (-6.8%, razon 0.50) | potencia 15948.6 (-7.1%, p = 1.7) | NIST 17162.5

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
