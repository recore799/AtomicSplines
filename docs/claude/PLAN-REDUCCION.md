# Plan de revisión y reducción del manuscrito

Escrito el 2026-09-24, a partir de la primera compilación completa (166 páginas). Objetivo: acercar
la tesis a unas 100 páginas y dejar el capítulo de resultados listo para un jurado de físicos. El
código sigue congelado (ver `CLAUDE.md`); todo lo de aquí es trabajo sobre el texto, salvo las
salidas nuevas de la etapa 3 que se indican.

**Punto de entrada de cada sesión mientras dure la reducción: este archivo.**

**Auditoría del 2026-10-09:** `REVISION-2026-10-09.md` tiene el estado del manuscrito de 82 páginas
y la lista de lo que hay que corregir antes de enviarlo, con archivo y línea. Las fases D, E y F
de abajo siguen pendientes y esa lista las concreta.

---

## 0. Dónde están las páginas

Páginas impresas de la compilación del 2026-09-24. El objetivo por bloque es orientativo.

| bloque | ahora | objetivo | palanca principal |
|---|---|---|---|
| preliminares (portada, índice, figuras) | 8 | 6 | menos subsecciones y figuras |
| prólogo + introducción + portadillas | 9 | 6 | fundir el prólogo en la introducción |
| cap. 2 partículas independientes | 13 | 9 | recortar prosa; corregir errores |
| cap. 3 correlación y estructura fina | 11 | 9 | notación y analogía de EDP del CI |
| cap. 4 B-splines y Galerkin | 9 | 5 | 1 de 3 figuras; STO/GTO aparece 3 veces |
| cap. 5 ciclo SCF | 13 | 9 | diagrama de flujo; tensores y caché al apéndice D |
| cap. 6 resultados | 28 | 18 | reordenar; §6.2; fusionar cuadros |
| cap. 7 discusión | 9 | 5 | repite el cap. 6 |
| cap. 8 conclusiones | 4 | 3 | repite §7.3 |
| apéndices A-C | 17 | 7 | unidades a un cuadro; Hund al cap. 3; 3j y CFP generales fuera |
| apéndice D | 19 | 4 | repite el cap. 4; dos listados de Julia; prosa de rendimiento |
| apéndice E | 17 | 4 | repite el cap. 4; prueba de Leibniz de 3 páginas |
| apéndice F | 7 | 5 | |
| bibliografía | 2 | 2 | |
| **total** | **166** | **~92** | |

**Formato adoptado el 2026-09-27** (la facultad no fija reglas; el usuario pidió uno estándar): A4,
12 pt, interlineado sencillo, márgenes de 2.5 cm y 3 cm del lado del lomo, encabezado en
minúsculas con `fancyhdr`, leyendas en `\small`, índice hasta secciones. Solo el formato bajó la
tesis de 167 a 129 páginas; con la fase B quedó en **116** (cuerpo 1-60, apéndices 61-109).
Desde esa fecha Claude compila y revisa el PDF (regla 1 de `CLAUDE.md`).

| bloque (tras la fase B; tras la C los apéndices bajaron a 18 y la tesis a 82) | páginas |
|---|---|
| preliminares | 5 |
| prólogo, introducción, portadillas | 7 |
| caps. 2-5 | 36 |
| cap. 6 resultados y discusión | 14 |
| cap. 7 conclusiones | 3 |
| apéndices A-F | 49 |
| bibliografía | 2 |

---

## Fase A. Correcciones de consistencia (no dependen de la estructura)

Afirmaciones que contradicen los cuadros, el código u otra parte del texto.

- [x] Intro: el 2 % "en los cuatro elementos" (solo Ge y Sn con α_d ajustado); "correlación estática".
- [x] Signo de los singletes: el ¹S₀ queda bajo en Si, Ge y Sn; la casi degeneración explica solo el
      carbono. Discusión (balance, sección de casi degeneración, déficit de mezcla), conclusiones.
- [x] Extrapolación de singletes: el ¹S₀ de Si converge peor y el ¹D₂ de Sn no admite extrapolación
      (razón 0.93); el texto decía que las dos estimaciones coincidían en todos.
- [x] C-DIIS: sin C-DIIS Ge y Sn sí convergen (241 y 250 iteraciones con level shift permanente de
      3 Ha); el texto decía que no convergían o que era impracticable.
- [x] V_pol: el apéndice F y resultados usaban W²/r⁴; el código y los caps. 3 y 5 usan W₆ sin
      elevar al cuadrado.
- [x] Frontera de Poisson en el cap. 4: quedaba Dirichlet; el código y los caps. 5 y E usan Robin.
- [x] Cap. 5: "repulsión apantallada" (V_pol es atractivo); "límite multiconfiguracional"; bandas de
      valencia y conducción.
- [x] Cuadratura: el código usa k + 6 puntos, no k; y solo S y T se integran sin error.
- [x] Costo lineal y matrices banda: el código almacena denso (Cholesky denso, sustituciones O(N²)).
- [x] Apéndice F: ARPACK y paralelización no están en el código (diagonalización densa de LAPACK).
- [x] Cap. 2: el intercambio es no local también en capa cerrada; "coordenada espaciotemporal".
- [x] Cap. 3: qué potencial entra en ξ(r) (el de Hartree de todo el átomo, más V_pol); ζ ∝ Z⁴ solo
      para iones hidrogenoides; números escritos a mano (∼30, ∼50 000, ≈1000 cm⁻¹); J = L + S;
      "electrodinámica cuántica"; la polarizabilidad no es un momento multipolar.
- [x] V_pol actúa también sobre los orbitales del core (entra en `ws.V`): documentado como limitación.
- [x] α_d = 0.75 en Ge frente a cruces en 0.66-0.67: la calibración elige sobre la rejilla.
- [x] El truncamiento en l_max = 3 no está acotado por las extrapolaciones en m (que son de m).
- [x] Menores: Si sin g del NIST, residuo de borrador en la discusión, "La Sección" por capítulo,
      etiqueta `tab:estanio_ci`.

Se dejan para la fase en la que se reescribe el pasaje, porque el pasaje se va o cambia entero:
§6.2 (α_d y la cola de valencia, el mínimo de V_rad), §6.8, apéndice A:72, apéndice B:49 (Z⁴),
apéndice C §C.5 (describe una evaluación de 3j con log-gamma y memoización que el código no tiene:
el CI usa `WignerSymbols.jl`), apéndice D "Estados espurios", apéndice E (O(N²) → O(N)).

## Fase B. Resultados y discusión (hecha el 2026-09-27)

Resultados y discusión quedaron en un solo capítulo, "Resultados y Discusión", con este orden:

1. Límite Hartree-Fock: energía y virial; F^k y ⟨r⁻³⟩ en un cuadro; Gram-Schmidt en un párrafo.
2. Espectro HF+CI de los cuatro elementos contra el NIST, con la convergencia en m resumida.
3. Del acoplamiento LS al intermedio: ζ frente a la separación (6/25)F² (no frente a F²), razón R,
   factor g y mezcla, con la estimación perturbativa ζ²/2Δ² que reproduce tanto la mezcla medida
   (con ζ_NIST y niveles del ASD) como la del modelo.
4. Déficit de ζ y V_pol en Ge y Sn, con la prueba cruzada de la ionización.
5. Limitaciones (antes §7.3), sin síntesis final propia: la hace el capítulo de conclusiones.

- [x] §6.2 a una figura y dos párrafos (fuera el pozo de V_rad y las figuras de potenciales,
      orbitales y densidades).
- [x] C-DIIS (cuadro y figura) al apéndice F; un párrafo en resultados.
- [x] Cuadros de convergencia (E_corr y singletes) al apéndice F; queda la figura de singletes.
- [x] Discusión fundida; la teoría de Zeeman y g, con la mezcla perturbativa, pasó al cap. 3
      (`sec:lande_teoria`). Fuera `discusion.tex`, §6.8, §7.4 y "Alcance del trabajo".
- [x] "Reproduce de manera cuantitativa una vez calibrada" eliminado.
- [x] Cuadros rediseñados en la etapa 3 (mismos números, comprobado contra la versión anterior):
      `parametros_radiales`, `ionizacion`, `niveles`, `niveles_vpol`, `barrido_vpol` (ζ en cm⁻¹),
      `acoplamiento` (nuevo), y encabezados agrupados en `convergencia_ci`, `cdiis`, `lande`.
      `scaling_law.pdf` grafica (6/25)F² en vez de F².
- [x] "Cuadro" → "Tabla" en el texto, como ya decían las leyendas.
- [ ] Figura de funciones de onda: tiene título interno y ejes poco legibles (etapa 4).

## Fase C. Apéndices (hecha el 2026-09-27: de 49 a 18 páginas; tesis en 82)

- [x] A: un párrafo y un cuadro de constantes y conversiones.
- [x] B (momento angular) eliminado: Hund y la regla de intervalos de Landé pasaron al cap. 3.
      Los apéndices se reletraron: A unidades, B tensores, C B-splines, D Galerkin, E auxiliares.
- [x] Tensores: fuera la evaluación de 3j con log-gamma (describía código que no existe) y los CFP
      generales; la fórmula del CI es ahora la general de parejas, la que implementa
      `np2_ci_full.jl` (normalización, fase de intercambio, 6j).
- [x] B-splines: nudos y condiciones de frontera, evaluación, separación geometría/estado (movida
      desde el cap. 5 §5.2, que conserva solo el operador de Poisson y sus fronteras) y el estado
      espurio, ahora medido (ver "Hallazgo" abajo). Fuera listados, BLAS y prosa de rendimiento.
- [x] Galerkin/Poisson: equivalencia Ritz-Galerkin (con la cota variacional), derivación directa de
      la ecuación de Poisson y condición de Robin. Fuera la formulación repetida, el cuadro de
      residuos ponderados y la prueba de Leibniz.
- [x] Auxiliares: Gram-Schmidt reescrito (sin "convergencia incondicional" ni Pauli), arquitectura
      en dos párrafos.

## Hallazgo del 2026-09-27: el espacio de orbitales del CI no es ortogonal al core

`build_orbital_pool` (np2_ci_full.jl) diagonaliza el Fock de core congelado sobre TODOS los
splines (1:n); el SCF usa 2:n-1 en s, 3:n-1 en p, 4:n-1 en d. En el canal s, el primer spline
(no nulo en r = 0) produce un autovalor espurio (-1178.6 Ha en C, -108515.4 en Sn) que el filtro
-Z^2 de `extract_virtuals` descarta, pero el resto del espectro s queda distorsionado y los
virtuales s solapan con los orbitales s del core hasta 3.6e-2 (peso máximo 2.0e-3 en un virtual).
p: <= 8.8e-4; d de Ge y Sn: <= 8.5e-3. Afecta solo a configuraciones con electrones s (¹S, ¹D), no
al ³P. Medido con `tesis/diagnostico_estados_espurios.jl`; los números salen en `valores_texto.md`.

Documentado como limitación (resultados §6.5, apéndice C). **No se corrigió**: cambia la física
del CI, que está congelada. Corregirlo = construir los virtuales sobre el mismo espacio que el SCF
(y/o ortogonalizarlos al core) y repetir la etapa 2 (~2.5 h, la corre el usuario). Decide el usuario.

## Fase D. Capítulos 2-5

- [ ] Prólogo: su párrafo útil (por qué construir desde cero) abre la introducción.
- [x] Intro: quitar la autoevaluación (":13, avala la capacidad técnica") y la de §8.4 (reescritas el 2026-10-09).
- [ ] Cap. 2: espín-estadística más corta; fuera la anécdota de la ocupación fraccionaria (2/3 F⁰).
- [x] Cap. 3: definir con precisión el modelo de CI (parejas nl n'l, m orbitales por canal, l_max).
- [x] Cap. 3: MP2 y la analogía de EDP fuera.
- [ ] Cap. 5: fuera el diagrama de flujo; §5.2 de tensores al apéndice D.

## Fase E. Estilo y notación

- [ ] "Conviene" (41 veces), "fundamental" (45) y el vocabulario inflado (21).
- [ ] V_eff tiene tres significados (cap. 4, cap. 5, resultados): unificar.
- [ ] Ecuaciones definidas en más de un lugar: expansión de Laplace (3), R^k (2), matrices de
      Galerkin (2), V_pol (3).

## Fase F. Procedencia de números (etapa 3)

Números del texto que no están en `valores_texto.md`. Agregarlos es de lo permitido sin preguntar.

- [x] Pisos del conmutador (`ciclo-scf.tex`, sección C-DIIS): en `valores_texto.md` desde el 2026-10-09,
      leídos de las trazas registradas.
- [ ] Prueba de malla del silicio (`resultados.tex`, 100/200/300 intervalos).
- [ ] Diferencias autovalor-Rayleigh (`discusion.tex`, limitaciones numéricas).
- [ ] Memoria del CI (apéndice F: 17 bytes por casilla, 1.1 y 1.7 GB).
- [x] ζ / (6/25)F² de HF por elemento: en el cuadro `acoplamiento`.
- [x] Peso perturbativo de ¹D₂ (ζ²/2Δ²) del modelo y del NIST: en `valores_texto.md`.

## Decisiones (resueltas el 2026-09-27)

- **Formato:** estándar elegido por Claude (ver §0); la facultad no fija reglas.
- **Fusionar cuadros:** autorizado; hecho en la fase B.
- **Estructura:** "Resultados y Discusión" en un solo capítulo; hecho.

Queda abierta la eliminación de las tres portadillas de parte (3 páginas).
