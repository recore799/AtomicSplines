# Plan de regeneración: del código a los números de la tesis

Escrito el 2026-09-13. Para el lado del código sustituye a la lista B1–B6 de `CIERRE-DEL-CODIGO.md`
y a las fases A y B de `FALTANTES-2026-09-13.md`. La redacción (fase C de `FALTANTES`) empieza
cuando termine lo de aquí.

---

## 0. En corto

- **Ningún `.jld2` registrado está mal.** Lo que se regenera es lo que sale de ellos: el CI, las
  tablas y las figuras. Además hacen falta 9 SCF nuevos con V_pol sobre orbitales ³P.
- **Tiempo de máquina:** ~30–40 min de SCF, ~2.5 h de CI y un par de minutos de tablas y figuras.
  Lo largo lo corres tú; el código se probó con corridas de humo baratas (sección 1).
- **Resultado:** 15 fragmentos de tabla en `docs/tablas/`, `valores_texto.md` con los números que
  cita la prosa y 11 figuras en `docs/figures/`, cada número con su procedencia.
- **Antes de correr, confirmar tres decisiones por defecto de `config.jl`:** orbitales ³P (2.1),
  energía con autovalores (2.2) y criterio de α_d (2.3).

---

## 1. Qué hay en el código

`benchmarks/np2_sequence/tesis/`:

| archivo | qué hace | probado |
|---|---|---|
| `config.jl` | único lugar de las decisiones y de los datos de literatura (NIST, Froese Fischer) | — |
| `comun.jl` | nombres de archivo, manifiesto, espacios de trabajo, `resultados_ci.toml`, formato de tablas | sí |
| `etapa1_scf.jl` | lista los `.jld2` que pide la configuración y corre los que faltan, con C-DIIS | lista: sí; las dos pruebas de `--prueba-vpol`: sí |
| `etapa2_ci.jl` | CI de valencia; guarda cada resultado al terminar y es reanudable | `--prueba` en C, Ge y Sn: sí |
| `etapa3_tablas.jl` | 15 fragmentos `tabular` y `valores_texto.md` | con el CI de prueba: sí |
| `etapa4_figuras.jl` | 11 figuras | con el CI de prueba: sí |
| `diagnostico_rayleigh.jl` | energía con autovalores frente a cociente de Rayleigh, sin SCF | corrido (2.2) |
| `analisis_zeta.jl` | ζ efectiva que pide el NIST, con y sin término tensorial | corrido (sección 6) |

Cambios en scripts existentes:

- **`germanium_rohf.jl` y `tin_rohf.jl` aceptan `alpha_d` y `r_c`.** Con `alpha_d = 0` el camino
  numérico es el de siempre. Con V_pol y C-DIIS reproducen los dos archivos registrados que
  generaron los scripts legados:
  - `germanium_rohf_results_av_R30.0_ad0.500.jld2`: a 1.35×10⁻⁸ Ha, en 72 iteraciones y 106 s;
  - `tin_rohf_results_3P_R30.0_ad1.000.jld2`: a 1.36×10⁻⁸ Ha, en 41 iteraciones y 246 s.
- **`germanium_rohf_vpol.jl` y `tin_rohf_vpol.jl` quedan como legado**, conservados como
  procedencia. El de germanio tenía invertido el signo de f₂ del ³P (−5/25, el mismo error de C y
  Si). Quedó corregido y nunca produjo un resultado registrado: su único archivo es `av`.
- **`verificar_resultados.jl`** lleva en su `CATALOGO` los 9 SCF planeados; `--emitir` los omite
  con un aviso mientras no existan.
- **La etapa 4 no pisa figuras a ciegas.** Los PDF no se versionan (`*.pdf` en `.gitignore`), así
  que la primera vez que reemplaza una figura copia la anterior a `docs/figures/anteriores/`.
- **Si un SCF no converge con C-DIIS en 300 iteraciones, la etapa 1 lo repite sola sin C-DIIS.**

---

## 2. Decisiones de `config.jl`, con los números que las sostienen

### 2.1 Orbitales ³P en todo (cambia la recomendación de `FALTANTES` §4.1)

`FALTANTES` recomendaba el promedio de configuración para el CI porque los barridos de V_pol
existentes son `av`. Al armar las tablas aparecieron dos argumentos en contra.

| | −ε_np ³P (eV) | −ε_np av (eV) | diferencia |
|---|---|---|---|
| C | 11.792 | 11.072 | 0.719 |
| Si | 8.085 | 7.585 | 0.500 |
| Ge | 7.819 | 7.329 | 0.491 |
| Sn | 7.212 | 6.764 | 0.448 |

- **El NIST mide la ionización desde el nivel fundamental ³P₀.** Con orbitales congelados, −ε de
  la Fock del ³P es la ionización desde el ³P; −ε de la Fock del promedio es la ionización desde
  el promedio de configuración. Difieren en ≈ (3/25) F², la diferencia entre los coeficientes de
  intercambio intra-capa (5/25 y 2/25).
- **El Cuadro `tab:koopmans` de la tesis (11.07 / 7.58 / 7.33 eV) usa el promedio.** Con ³P el
  carbono pasa a 11.79 eV y sobrestima el NIST (11.26).
- **Mezclar convenciones se filtra a cinco tablas y dos figuras.** Con una sola, las tablas que
  validan contra Froese Fischer (energía, ⟨r⁻³⟩, F^k), el CI y V_pol hablan de los mismos orbitales.
- **En ζ la elección pesa 1.5–1.8 %:** ζ(³P)/ζ(av) = 1.018, 1.015, 1.016 y 1.015.
- **Costo:** 9 SCF (~30–40 min). Con `ESTADO = "av"` y `ESTADO_SENSIBILIDAD = "3P"` bastan 3
  (germanio con α_d = 0.25, 0.75 y 1.0; ~8 min), porque el barrido del estaño en promedio ya existe.

### 2.2 Energía: autovalores por defecto (decisión tuya)

`diagnostico_rayleigh.jl` reconstruye la Fock de cada canal con los orbitales guardados y evalúa
las dos fórmulas sin correr ningún SCF. Sus tres controles (energía guardada, orbitales de core,
eigenvalor de valencia) dan ≤ 1×10⁻⁹ Ha.

| archivo | E_Rayleigh − E_autovalor (Ha) |
|---|---|
| C ³P | 5.3×10⁻¹² |
| Si ³P | −1.72×10⁻⁶ |
| Ge ³P | −4.3×10⁻⁷ |
| Sn ³P | −5.8×10⁻⁷ |

En promedio de configuración llega a −4.2×10⁻⁶ Ha (silicio).

- **A las cifras de Froese Fischer no se ve, pero la tabla imprime 8 decimales.** Con Rayleigh, el
  silicio pasaría de −288.85435715 a −288.85435887 y su Δ de 2.8×10⁻⁶ a 1.1×10⁻⁶; germanio y
  estaño se mueven en la séptima cifra decimal.
- **Por qué autovalores por defecto:** es el `E_total` que guardan los `.jld2` y registra
  `RESULTS.toml`, de modo que la tabla y el manifiesto dicen el mismo número, y es el que ya
  imprime la tesis.
- **Por qué Rayleigh es defendible:** es la energía de los orbitales guardados. La fórmula con
  autovalores deja de serlo en cuanto Gram-Schmidt saca al orbital de valencia de su autovector.
- **Cambiarlo es una línea** (`ENERGIA_HF = :rayleigh`) y no obliga a correr ningún SCF: el punto
  fijo lo definen las Fock, no la fórmula de la energía. Lo que se elija se declara en el Apéndice F.
- **La cifra de 1.29×10⁻⁵ Ha de `HALLAZGOS` §6 no aplica:** se midió con el signo de f₂ todavía
  invertido.

### 2.3 α_d: barridos y criterio automático

- **Barridos:** germanio α_d ∈ {0.25, 0.5, 0.75, 1.0} (el 0.5 de la tesis quedaba en el borde de
  su barrido viejo); estaño α_d ∈ {1, …, 6}. El radio de corte es r_c = 1.0 a₀ en los dos.
- **Criterio:** la columna CI+V_pol usa el α_d del barrido que minimiza el error relativo rms de
  ³P₁ y ³P₂ (`CRITERIO_ALPHA = :rms_3P`). El «3.5 de compromiso» de la tesis desaparece (no hay
  SCF con ese valor): se reporta el mejor punto y la tensión entre ³P₁ y ³P₂ se discute con la
  tabla completa del barrido.
- **Si el elegido cae en un extremo, el barrido no lo encierra:** hay que ampliarlo en `config.jl`
  y volver a la etapa 1.
- **Referencia gruesa, solo con ζ y sobre orbitales del promedio:** el ζ que pide el NIST se
  alcanza hacia α_d ≈ 0.8 en germanio (extrapolado más allá de su único punto, 0.5) y ≈ 3.3 en
  estaño. Con ³P, un poco menos.

### 2.4 Espacio activo y formato

- **Espacio activo:** m = 4, 8, 12, 16 y 20 orbitales por cada l ≤ 3. El último punto de la curva
  de convergencia es el HF+CI de las tablas, el mismo espacio del log del 07-09.
- **Formato:** fragmentos `tabular` sin leyenda ni etiqueta; esas se quedan en el capítulo.

---

## 3. Qué correr en tu máquina, en orden

Todo desde la raíz del repo.

| # | comando | tiempo | revisar |
|---|---|---|---|
| 0 | `julia benchmarks/verificar_resultados.jl` | s | 16/16 |
| 1 | `julia --project=. benchmarks/np2_sequence/tesis/etapa1_scf.jl` | s | 9 en `falta` |
| 1b | (opcional) `... etapa1_scf.jl --prueba-vpol` | ~6 min | los dos d ≈ 1.4×10⁻⁸ Ha; ya lo corrí el 13-09 |
| 2 | `julia --project=. benchmarks/np2_sequence/tesis/etapa1_scf.jl --correr 2>&1 \| tee etapa1.log` | ~30–40 min | cada SCF termina en `-> E = ...` |
| 3 | `julia --project=. benchmarks/verificar_resultados.jl --emitir > benchmarks/RESULTS.toml` y luego `julia benchmarks/verificar_resultados.jl` | s | 25/25; el diff trae solo 9 entradas nuevas |
| 4 | commit de los 9 `.jld2` y el manifiesto | — | — |
| 5 | `julia --project=. benchmarks/np2_sequence/tesis/etapa2_ci.jl --lista` | s | ninguno `SIN REGISTRAR` |
| 6 | `julia --project=. benchmarks/np2_sequence/tesis/etapa2_ci.jl 2>&1 \| tee etapa2.log` | ~2.5 h | cada archivo pasa T1–T4 (`TODO OK`) |
| 7 | commit de `benchmarks/np2_sequence/resultados_ci.toml` | — | — |
| 8 | `julia --project=. benchmarks/np2_sequence/tesis/etapa3_tablas.jl` | ~15 s | termina en «ninguna celda quedó sin dato» |
| 9 | `julia --project=. benchmarks/np2_sequence/tesis/etapa4_figuras.jl` | ~1 min | 11 figuras; las anteriores quedan en `docs/figures/anteriores/` |
| 10 | commit de `docs/tablas/` | — | las figuras no se versionan |

Notas:

- **Paso 6 en tandas:** `--grupo convergencia` (~25 min), `sensibilidad` (~20 min) y `barrido`
  (~90 min), o por elemento con `--solo`. Si se interrumpe, relanzar el mismo comando: lo ya
  guardado no se repite.
- **`resultados_ci.toml` es la procedencia:** los logs (`*.log`) están ignorados por git y no hacen
  falta para reconstruir nada.

---

## 4. Lo primero que hay que revisar cuando lleguen los números

1. **Reproducibilidad.** A m = 20, el HF+CI sobre los `_3P_` sin V_pol tiene que reproducir el log
   del 07-09, que usó los mismos archivos y el mismo espacio:
   - C: E_corr = −8.698045×10⁻³ Ha, ³P₂ = 63.08 cm⁻¹
   - Si: −8.562860×10⁻³ Ha, 210.58 cm⁻¹
   - Ge: −7.060409×10⁻³ Ha, 1257.58 cm⁻¹
   - Sn: −6.244565×10⁻³ Ha, 2921.70 cm⁻¹

   Si no coinciden, parar e investigar antes de escribir nada.
2. **Autoverificación y raíces.** T1–T4 en `etapa2.log` para cada archivo, y `raices = [1, 1, 1]`
   en `resultados_ci.toml`. Una raíz distinta de 1 no es un error: es información física que va a
   la tesis.
3. **Convergencia.** Carbono es el único con cola apreciable (~1.5×10⁻³ Ha según `HALLAZGOS`).
4. **Sensibilidad ³P frente a promedio.** Los ³P_J deben diferir ~1.6 %, la razón de ζ. Cuánto se
   mueven ¹D₂ y ¹S₀ no se ha medido, y es lo que dice si la elección de orbitales importa para los
   singletes.
5. **Barridos.** Que el α_d elegido no esté en un extremo, y cuánta tensión queda entre ³P₁ y ³P₂
   en estaño.
6. **Koopmans.** Con orbitales ³P el carbono sobrestima el potencial de ionización (−ε = 11.79 eV
   frente a 11.26 del NIST), y la corrección CI lo empuja más arriba. Koopmans no incluye la
   relajación del ion, que baja el potencial, y sumar solo la correlación del neutro rompe esa
   cancelación parcial.

---

## 5. De los números al manuscrito

Cada tabla del capítulo se reemplaza por su fragmento, conservando leyenda y etiqueta:

```latex
\begin{table}[htbp]
    \centering
    \input{tablas/energia_global}
    \caption{...}
    \label{tab:resultados_energia_global}
\end{table}
```

| etiqueta actual | fragmento |
|---|---|
| `tab:resultados_energia_global` | `energia_global.tex` |
| `tab:momentos_inversos` | `momentos_inversos.tex` |
| `tab:integrales_slater` | `integrales_slater.tex` |
| `tab:koopmans` | `koopmans.tex` |
| `tab:niveles_energia_ligeros` | `niveles_ligeros.tex` |
| `tab:niveles_energia_pesados` | `niveles_pesados.tex` |
| `tab:estanio_ci` | `barrido_vpol_Sn.tex` (y el nuevo `barrido_vpol_Ge.tex`) |
| `tab:germanio_ionizacion` | `ionizacion_Ge.tex` (hay también `ionizacion_Sn.tex`) |
| nuevas | `convergencia_ci.tex`, `convergencia_ci_costo.tex`, `cdiis.tex`, `zeta.tex`, `lande.tex` |

Figuras:

- **Las siete actuales se regeneran con el mismo nombre.**
- **Nuevas:** `convergencia_ci.pdf`, `cdiis_convergencia.pdf`, `zeta_razon.pdf` y `barrido_vpol.pdf`.

Los números que la prosa cita fuera de las tablas están en `docs/tablas/valores_texto.md`:

- mínimos de V_rad, incluido el del estaño que faltaba (−399.43 Ha en r = 0.042 a₀);
- ζ y F² en cm⁻¹, y los ajustes de ζ al NIST;
- E_corr, g_eff y niveles;
- el error de cada punto de los barridos.

Hallazgos de esta regeneración que cambian texto o tablas (además de `FALTANTES` §3.2):

- **El virial del estaño que imprime la tesis (2.00000000) no sale del archivo registrado.** Con
  la misma fórmula de `tin_validate.jl`, y con la E y la T que la propia tabla ya imprime, da
  1.99999938. El renglón del germanio también difiere en la última cifra de T (2075.36256 impreso,
  2075.36257 calculado).
- **Dos figuras publicadas tenían la base mal construida.** `np2_all_orbitals.pdf` y
  `np2_total_densities.pdf` salían de `examples/scratch/plot_all_*.jl`, que usaba
  `N_elems = num_splines − K` (un elemento de menos). Evaluaba los coeficientes sobre una malla de
  nudos que no era la del cálculo.
- **`np2_radial_potentials.pdf` cortaba al estaño:** el eje fijo en −300 Ha deja fuera su pozo.

---

## 6. El problema de ζ: qué es y si hace falta otro modelo

Valores de esta sesión (ζ en cm⁻¹, orbitales ³P, ajustes de `analisis_zeta.jl`). La etapa 3 los
regenera en `zeta.tex` y `valores_texto.md`.

| | (αZ)² | ζ núcleo desnudo | ζ calculado | ζ NIST (a) | ζ NIST (b) | calc / NIST (a) | R NIST | R modelo |
|---|---|---|---|---|---|---|---|---|
| C | 0.0019 | 59.3 | 42.2 | 28.9 | 32.6 | 1.460 | 1.646 | 1.983 |
| Si | 0.0104 | 167.5 | 141.0 | 148.1 | 148.5 | 0.952 | 1.895 | 1.911 |
| Ge | 0.0545 | 896.3 | 823.7 | 920.6 | 910.4 | 0.895 | 1.531 | 1.563 |
| Sn | 0.1331 | 1992.5 | 1885.9 | 2241.1 | 2207.5 | 0.842 | 1.026 | 1.137 |

Cómo se leen las columnas:

- **ζ núcleo desnudo** = (α²/2) Z ⟨r⁻³⟩ con el mismo orbital: el ζ de un electrón que viera la carga
  nuclear completa.
- **ζ calculado** = `compute_zeta`: (α²/2) ⟨(1/r) dV/dr⟩, con V = −Z/r más el potencial de Hartree
  de todos los electrones, incluido el propio (`compute_effective_central_potential`).
- **ζ NIST.** (a) es la ζ única que, con energías LS libres, reproduce ³P₂, ¹D₂ y ¹S₀ en las
  matrices de Breit-Pauli de p²; (b) reproduce ³P₁, ¹D₂ y ¹S₀.
- **R** = (³P₂ − ³P₁)/(³P₁ − ³P₀). Una ζ de un cuerpo en LS puro da R = 2 (regla de Landé), y el
  acoplamiento intermedio lo baja. «R modelo» sale del CI a m = 20 del 07-09.
- **`HALLAZGOS` §5 daba 917.1 y 2168.8** para Ge y Sn con otra extracción; la tesis tiene que decir
  cuál usa.

### ¿Está mal el Breit-Pauli? No

1. **Límite jj.** Las matrices reproducen el espectro jj cuando ζ domina; la tesis ya lo muestra
   en `correlacion-estructura-fina.tex:120`.
2. **El propio experimento valida la estructura.** Con una sola ζ ajustada a tres niveles del
   NIST, el cuarto se predice a 0.3 % en silicio (³P₁ = 76.9 frente a 77.1), a 1.3 % en germanio
   (564.5 frente a 557.1) y a 2.1 % en estaño (1727.3 frente a 1691.8). La forma del modelo es la
   correcta; lo que falla es el valor de ζ.
3. **La excepción es el carbono** (14.5 frente a 16.4, −12 %), y ahí el experimento muestra física
   que ninguna ζ de un cuerpo contiene. El NIST da R = 1.646 donde el modelo da 1.983:
   - dentro de un término LS, todo operador de rango (1,1) (espín-órbita y espín-otra-órbita)
     cumple la regla de Landé, y solo el espín-espín, de rango (2,2), la rompe;
   - con un término tensorial en el ajuste, los cuatro niveles salen exactos con ζ_ef = 28.0 y
     D = −0.31 cm⁻¹.

### De dónde sale la discrepancia

- **Ligeros (C y Si): el ζ calculado sobra.** `compute_zeta` apantalla el término de un cuerpo con
  el potencial de Hartree, pero no contiene el intercambio ni el espín-otra-órbita de los términos
  de dos cuerpos de Breit-Pauli (el ζ de Blume-Watson). En la práctica esos términos lo reducen
  más, y pesan como ~1/Z_ef frente al término nuclear. El ζ(2p) = 31.946 cm⁻¹ de Froese Fischer que
  usaba `carbon_fine_structure.jl` cae entre el nuestro (42.2) y el del NIST (28.9): si es el de
  Blume-Watson (verificarlo en el libro), es la prueba directa.
- **Pesados (Ge y Sn): el ζ calculado falta, y ningún apantallamiento lo arregla.**
  - Para un potencial de Hartree, dV/dr = (Z − q_enc(r))/r² ≤ Z/r², así que ζ calculado ≤ ζ
    desnudo con el mismo orbital; la tabla lo cumple en los cuatro.
  - El NIST pide 2.7 % (Ge) y 12.5 % (Sn) *más* que esa cota. Con estos orbitales no relativistas,
    ninguna forma de apantallar este ζ llega al experimento: hace falta un orbital de valencia más
    contraído cerca del núcleo.
  - Hay dos mecanismos conocidos que lo contraen, y con dos elementos los datos no los separan: la
    **contracción relativista**, de orden (αZ)² = 5 % y 13 %, y la **polarización del core**, que
    es justo lo que modela V_pol (sube ζ de 810.7 a 882.4 cm⁻¹ en germanio con α_d = 0.5, sobre
    promedio de configuración).
- **El cruce monótono** de calc/NIST con Z (1.46 → 0.84) es la suma de las dos omisiones de signo
  contrario, cada una dominante en un extremo.

### ¿Hace falta cambiar el modelo para la tesis? No

El diagnóstico se demuestra con lo que el código ya produce (`zeta.tex`, `zeta_razon.pdf` y
`valores_texto.md`), y cada mejora es un proyecto propio:

| mejora | arregla | costo |
|---|---|---|
| ζ de Blume-Watson (intercambio y espín-otra-órbita con orbitales HF) | C y Si | integrales radiales nuevas sobre la maquinaria de Poisson; días, más validación |
| término espín-espín | R del carbono | similar al anterior |
| orbitales relativistas: Dirac-Hartree-Fock con B-splines, o las correcciones escalares en la ecuación radial del HFR de Cowan (ya citado en la tesis) | Ge y Sn | otro código; semanas |

**Lo que sí hay que escribir.** Calibrado contra el NIST, α_d absorbe también la parte del déficit
que no es polarización (la relativista), así que en Ge y Sn es un parámetro efectivo y no la
polarizabilidad del core. Compararlo con la polarizabilidad dipolar de la literatura para esos
cores dice cuánto de él es polarización real.

---

## 7. Legado

Estos archivos no se borran: son procedencia o historia. Ninguno alimenta ya un número de la tesis.

- **Generadores viejos de tablas y figuras:** `generate_table.jl`, `calibrate_vpol.jl`,
  `calibrate_tin_vpol.jl`, `scaling_law.jl` (constantes a mano), `plot_energy_levels.jl` (niveles a
  mano) y `examples/scratch/plot_*.jl` (archivos sin registrar, base equivocada).
- **Scripts de estructura fina previos al CI:** `*_fine_structure.jl` y `lande_g_factors.jl`.
- **Scripts de V_pol:** `*_rohf_vpol.jl`, conservados como procedencia.
- **El CI de pares de `np2_toy_ci.jl`.** Sus utilidades (`compute_zeta`, `get_h_core`) las sigue
  usando `np2_ci_full.jl`.

**Tabla de C-DIIS.** Sale de las trazas existentes, con el número de iteraciones. Las trazas sin
C-DIIS de Ge y Sn se reconstruyeron del log impreso, con E a 8 decimales: alcanza para contar
iteraciones, y la tabla marca como «< 10⁻⁸» la diferencia de energía que queda por debajo de esa
precisión. Si la tesis quiere tiempos de pared o más cifras, `diis_benchmark.jl` (~45 min) vuelve
a escribir las ocho trazas.

---

## 8. Cuándo se congela el código

- **Manifiesto:** `verificar_resultados.jl` da 25/25.
- **CI:** `etapa2_ci.jl --lista` marca todo como `listo`.
- **Tablas:** `etapa3_tablas.jl` termina sin faltantes. No cuentan las celdas `--` que dependen de
  literatura: F² de referencia del Ge, virial de referencia del Sn, g del NIST de C, Si y Sn.
- **Revisión:** la sección 4 está hecha.
- **Commit:** manifiesto, `resultados_ci.toml` y `docs/tablas/`.

A partir de ahí no se toca código. Una discrepancia que aparezca al redactar se documenta como
limitación.
