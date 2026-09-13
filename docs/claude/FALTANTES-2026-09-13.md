# Qué falta para terminar la tesis — inventario del 2026-09-13

Contrasta `CIERRE-DEL-CODIGO.md` (08-09) y el resumen del estado de datos del 07-09 contra el
repo tal como está hoy. La verificación del LaTeX fue estática, sin compilar. Se corrió
`verificar_resultados.jl` y un diagnóstico de solo lectura en Julia: ζ y ⟨r⁻³⟩ de los 16
resultados registrados, con el `compute_zeta` de `np2_toy_ci.jl`, que reproduce los ζ del log
del 07-09 a menos de 2×10⁻⁹ Ha. No se tocó ningún número de la tesis ni ningún `.jld2`.

> **Actualización del mismo día.** Para el lado del código (secciones 2 y 4, y las fases A y B de
> la sección 5) manda `PLAN-REGENERACION.md`, que además cambia la recomendación de §4.1: ahora
> los orbitales son ³P en todo. La sección 3, el manuscrito, sigue vigente.

---

## 0. En corto

- **El código no está congelado.** De los seis bloqueantes de `CIERRE`, solo B2 está cerrado.
  B1 no ha empezado, y B3 tal como está definido no alcanza: `generate_table.jl` cubre dos de
  las ocho tablas del capítulo de resultados (y las funde en una) y ninguna figura.
- **La física de B1 se confirma con números** (V_pol sube ζ), pero la narrativa invertida
  también está en el capítulo de teoría, no solo en `resultados.tex`.
- **Dos figuras publicadas salen de constantes escritas a mano y cinco de archivos fuera del
  manifiesto.**
- **Falta escribir:** conclusiones (vacío), trabajo futuro (no existe), C-DIIS (cero líneas) y
  la evidencia de convergencia del CI (la curva no se guardó). **Falta reescribir:** la
  discusión completa, `sec:ci_method` y la introducción.
- **Hay decisiones que van antes de regenerar** (sección 4). La que más pesa: qué orbitales,
  promedio o ³P, alimentan el CI y Breit-Pauli.
- Hay dos commits sin subir.

---

## 1. Estado del repositorio

- **Push pendiente.** `origin/main` está en `352b4ed`; faltan `63b2555` y `31e5e2c`, y este
  último trae `tin_rohf_results_av_R30.0.jld2`.
- **Sin trackear:** `docs/claude/CIERRE-DEL-CODIGO.md`.
- **El log de la tabla está ignorado.** `docs/claude/tabla_np2_2026-09-07.log` cae en `*.log`
  (`.gitignore:30`) y cita los nombres `_f2fix`, que ya no existen, así que no puede ser la
  procedencia final. Sus números sí sirven como vista previa para redactar: los ζ de los
  archivos canónicos los reproducen. Lo que no sirven es para pegarse en las tablas.
- **`docs/chapters/#resultados.tex#` no es trabajo perdido.** Es un autosave de Emacs del 05-09
  a las 18:20. Sus 53 párrafos que no aparecen en el manuscrito actual son versiones anteriores
  del mismo texto (registro inflado, tres elementos en lugar de cuatro). Se puede borrar.
- **No están en el repo** los documentos hermanos que cita el resumen del 07-09
  (`estado_codigo_ci.md`, `estado_manuscrito.md`) ni ese mismo resumen.
- **Manifiesto:** `verificar_resultados.jl` da 16/16 hoy.
- **Datos fuera del manifiesto:** cinco figuras de resultados leen
  `examples/scratch/{carbon,silicon,germanium}_rohf_results_R30.0.jld2`, del 24-jun, sin
  trackear y sin sufijo de estado. Es el mismo riesgo que la reorganización eliminó en
  `benchmarks/`.

---

## 2. Bloqueantes de `CIERRE-DEL-CODIGO.md`: estado real

| | Estado | Evidencia |
|---|---|---|
| B1 | abierto, sin empezar | `calibrate_vpol.jl:44` y `calibrate_tin_vpol.jl:47` siguen llamando a `run_toy_ci` (legado); `file_vpol = nothing` en los cuatro casos de `generate_table.jl` |
| B2 | **cerrado** | commits `6e19da8` y `4716d4c`; no queda ninguna referencia a `_f2fix` |
| B3 | abierto | la corrida del 07-09 usó los `_f2fix`, sin V_pol, y su log está ignorado; además el alcance no basta (2.2) |
| B4 | a medias | el archivo existe (sin push), pero la columna sigue en `{-}` (`resultados.tex:241-245`) y falta la decisión 4.1 |
| B5 | abierto, más grande de lo que dice | 2.3 |
| B6 | abierto | 2.4 |

### 2.1 B1: lo que confirma el diagnóstico

| archivo | α_d | ζ (cm⁻¹) | ⟨r⁻³⟩ (a₀⁻³) |
|---|---|---|---|
| C av | 0 | 41.480 | 1.661837 |
| C ³P | 0 | 42.216 | 1.691811 |
| Si av | 0 | 138.903 | 2.016763 |
| Si ³P | 0 | 140.992 | 2.047231 |
| Ge av | 0 | 810.736 | 4.717950 |
| Ge ³P | 0 | 823.671 | 4.793275 |
| Ge av | 0.5 | 882.355 | 5.130609 |
| Sn av | 0 | 1858.739 | 6.720989 |
| Sn ³P | 0 | 1885.936 | 6.819372 |
| Sn ³P | 1 | 1978.220 | 7.152685 |
| Sn av | 1 | 1950.160 | 7.051188 |
| Sn av | 2 | 2062.431 | 7.456421 |
| Sn av | 3 | 2196.675 | 7.940721 |
| Sn av | 4 | 2354.268 | 8.509012 |
| Sn av | 5 | 2536.771 | 9.166898 |
| Sn av | 6 | 2745.869 | 9.920414 |

El ⟨r⁻³⟩ se integró por trapecio sobre `R_grid`: es diagnóstico, no número de tabla.

1. **V_pol sube ζ** cuando se compara dentro de una misma convención. En germanio (promedio),
   ζ pasa de 810.7 a 882.4 cm⁻¹ con α_d = 0.5 (+8.8 %); en estaño (promedio), de 1858.7 a 2196.7
   con α_d = 3 (+18.2 %). El orbital se contrae (⟨r⁻³⟩ de 4.718 a 5.131 en Ge), como dice
   `CIERRE`. El texto de la tesis dice lo contrario.
2. **La elección entre ³P y promedio mueve ζ entre 1.5 y 1.8 %**, no el ~4 % de `HALLAZGOS` §5:
   la razón ζ(³P)/ζ(av) vale 1.018 (C), 1.015 (Si), 1.016 (Ge) y 1.015 (Sn). El Cuadro
   `tab:niveles_energia_pesados` mezcla ³P (HF+CI) con promedio (CI+V_pol). Eso *subestima* el
   efecto de V_pol en unos 1.7 puntos; no lo fabrica.
3. **Estimación por escalamiento lineal con ζ** (no es número de tabla): el ³P₂ del germanio
   quedaría en ≈ 1238 cm⁻¹ con HF+CI sobre promedio (−12 %) y en ≈ 1347 con CI+V_pol (−4.5 %),
   frente a 1410 del NIST. La conclusión de `CIERRE` sobrevive.
4. **El CI sí ve V_pol.** Los `_ad*.jld2` guardan `alpha_d` y `r_c`, y `run_full_ci` los pasa a
   `cached_init_scf_workspace` (`np2_ci_full.jl:389-415`). No hay error silencioso.
5. **La parte de CI de B1 no necesita SCF nuevos.** Los SCF de Ge con α_d = 0.5 y de Sn con
   α_d = 1…6 existen y están registrados. Solo hacen falta SCF para un α_d fuera del barrido (el
   3.5 de `resultados.tex:292` no tiene archivo) o para un barrido de germanio.
6. **El alcance textual de B1 es mayor** que `resultados.tex:248, 263-265, 292`:
   - `correlacion-estructura-fina.tex:163-167` sostiene el mecanismo equivocado (relajación
     hacia afuera, ⟨1/r³⟩ menor, HF que sobrestima ζ) y contradice a `:154` del mismo capítulo,
     que dice lo correcto: V_pol empuja al electrón de valencia hacia el núcleo.
   - `ciclo-scf.tex:192` habla de «repulsión apantallada», pero V_pol es atractivo.
   - En `resultados.tex:269`, V_pol casi no toca los singletes; sin embargo, en los números
     viejos el ¹S₀ se movía +699 cm⁻¹. En `:271`, el error del ¹S₀ es «superior a los
     3000 cm⁻¹»; con el CI nuevo, HF+CI da 15780 frente a 16367 (−587). Además, `:271` atribuye
     el residuo a «ampliar el espacio de configuraciones», lo que contradice la decisión 1 de
     `CIERRE`.
   - `resultados.tex:248` da α_d = 0.5 para el germanio sin decir de dónde sale. Si se calibró,
     hay que decirlo y mostrar el barrido; si viene de la literatura, falta la cita.
   - El comentario de `np2_toy_ci.jl:415` repite el mecanismo al revés (legado, cosmético).

### 2.2 B3: por qué `generate_table.jl` no alcanza

La regla 2 de `CLAUDE.md` y el propio `CIERRE` («ningún número del capítulo 6 debe sobrevivir
sin pasar por aquí») piden procedencia para todo el capítulo, y hoy casi nada la tiene.

| objeto | generador actual | problema |
|---|---|---|
| `tab:resultados_energia_global` | ninguno (solo `tin_validate.jl` calcula el virial de Sn) | renglón de Sn mal formado (`:40-42`: la referencia cae en la columna Δ y faltan las Δ); `:11` no menciona Sn |
| `tab:momentos_inversos` | ninguno (solo Sn, en `tin_validate.jl`) | reproduce desde los `_3P_` registrados a ≤10⁻⁵; solo falta la procedencia |
| `tab:integrales_slater` | ninguno (el T3 del CI imprime F², no F⁰) | B6 |
| `tab:koopmans` | ninguno | columna CI contaminada (2.3); renglón de Sn vacío |
| `tab:niveles_energia_ligeros` y `_pesados` | `generate_table.jl` | emite una sola tabla de cuatro elementos (`tab:niveles_energia`), no las dos del manuscrito; `file_vpol` vacío |
| `tab:estanio_ci` | `calibrate_tin_vpol.jl` (CI legado) | la columna ζ reproduce (0.00940…0.01251); ³P₁ y ³P₂ salen del CI legado |
| `tab:germanio_ionizacion` | ninguno | B5 |
| `fig:scaling_law` | `scaling_law.jl` | constantes a mano (`:23-24`): ζ(C) = 31.9, ζ(Ge) = 930, ζ(Sn) = 2276; F²(Si) = F²(Ge) = 0.1627 |
| `fig:np2_levels` | `plot_energy_levels.jl` | niveles a mano (`:142-145`), los del CI viejo; el estaño lleva números que no están en ninguna tabla ni tienen procedencia; el caption habla solo del germanio |
| `fig:wavefunctions`, `fig:all_orbitals`, `fig:total_densities`, `fig:potentials_eff`, `fig:potentials_rad` | `examples/scratch/plot_*.jl` | leen los archivos del 24-jun fuera del manifiesto para C, Si y Ge y el `_3P_` para Sn (convención mezclada); de ahí salen los mínimos de V_rad de `:97`, y falta el del estaño (`\todo` de `:98`) |

Texto que depende de esos números:
- `:97-98`.
- `:172-174`, que usa el ζ de Froese Fischer: «~30 cm⁻¹» para C frente a 42.2 calculado, y
  «cerca de 1000» para Ge frente a 824.
- `:182`, donde el F² de Si pasa de 0.163 a 0.166.

**B3 tiene que ser una pasada única que emita las ocho tablas y las siete figuras desde
archivos registrados.** Rehacer solo `generate_table.jl` no basta.

### 2.3 B5: también la tabla de Koopmans

La columna «CI (ΔE)» de `tab:koopmans` (`resultados.tex:198-200`) es, número por número,
−ε_np + |E_corr| del CI legado: +0.021 eV (C), +0.128 eV (Si) y +0.075 eV (Ge). Son los
−7.75×10⁻⁴, −4.72×10⁻³ y −2.76×10⁻³ Ha de `PLAN` §1. El problema no se limita al renglón de
germanio de `tab:germanio_ionizacion`.

Con el E_corr del CI nuevo la corrección sería de +0.24, +0.23 y +0.19 eV. Como estimación, no
número de tabla, eso da C ≈ 11.31 eV, Si ≈ 7.81 eV y Ge ≈ 7.52 eV, frente a 11.26, 8.15 y 7.90
del NIST. De ahí se siguen varias cosas:
- `:208` («apenas centésimas de electrón-volt») deja de ser cierto.
- El carbono pasaría a sobrestimar, y hay que explicarlo. Koopmans no incluye la relajación del
  ion, que baja el potencial de ionización; sumar solo la correlación del neutro rompe esa
  cancelación parcial.
- Cambian las tres filas de `tab:germanio_ionizacion`, porque la de V_pol también usa el CI.
- `:311` se reescribe con la decisión 1 de `CIERRE`: el espacio está convergido y el residuo es
  del modelo.
- Falta el renglón de Sn en `tab:koopmans`.

### 2.4 B6: integrales de Slater

- **Renglón de Si:** van los valores recalculados, 0.32969971 y 0.16594819 (`HALLAZGOS` §2); el
  F² lo confirma el T3 del log del 07-09.
- **F⁰ de Ge:** el 0.30137782 impreso no reproduce con ningún cálculo; el recalculado es
  0.31650266.
- **Referencia F.F. de Ge:** el F²(4p,4p) es una copia del de silicio (0.16593265). Hay que
  recuperarlo del libro de Froese Fischer (1977); el código no puede darlo.
- **Sn:** los valores impresos están truncados (0.27874800 y 0.14723200 frente a 0.27874830 y
  0.14723221), así que su Δ cambia en la tercera cifra.
- **Texto:** `resultados.tex:178` atribuye el error de Si a la discretización y
  `discusion.tex:14` al promedio de configuración. Ambas atribuciones son falsas: la causa fue
  el signo de f₂ más una transcripción equivocada.

---

## 3. Manuscrito

### 3.1 Falta por completo

- **Conclusiones.** `conclusiones.tex` es un marcador de posición de tres palabras.
- **Trabajo futuro.** No existe en ningún lado, y `CIERRE` le asigna contenido:
  - el cruce de ζ_calc/ζ_NIST con Z como resultado en sí mismo;
  - un CI multirreferencia con el ns activo para el ¹S;
  - los términos de dos cuerpos de Breit;
  - la contracción relativista (Dirac-Hartree-Fock con B-splines; `Johnson1988` ya está en la
    bibliografía).

  Recomendación: una sección dentro de Conclusiones.
- **C-DIIS.** Cero líneas en el manuscrito, aunque los datos existen (`HALLAZGOS` §6 y los
  `diis_trace_*.csv`, trackeados). Según `PLAN` §6:
  - la formulación va en `ciclo-scf.tex`, después del level shifting (`:41-59`);
  - la tabla de iteraciones y la curva del residual van en `resultados.tex`;
  - la implementación va en el Apéndice F: historial, cond(B), cambio de régimen por residual,
    salvaguardas, y el piso del conmutador causado por los multiplicadores de Gram-Schmidt,
    con su proyección.
- **Convergencia del CI.** `PLAN` (fase 1) dice que «va a la tesis»: la curva de E_corr contra el
  tamaño del espacio, con orbitales, CSFs, R^k y tiempo de pared. **La curva no se guardó.**
  `HALLAZGOS` conserva solo el renglón m = 20 y las colas, y `examples/scratch/ci_convergence.csv`
  es del 24-jun, con el CI legado. Hay que volver a correr `convergence_study` guardando un CSV:
  del orden de media hora para los cuatro elementos.
- **Límites del modelo, declarados explícitamente:**
  - la casi-degeneración ns² ↔ np², que deprime al ¹S (decisión 4 de `CIERRE`);
  - la cancelación core-core en los desdoblamientos (decisión 2);
  - el término de dos cuerpos de V_pol que se omitió (Apéndice F, `:112`), que hay que conectar
    con la discusión de ¹D y ¹S.
- **Preliminares:** no hay resumen ni agradecimientos, y `\maketitle` es genérico. Revisar el
  formato que pide la institución.

### 3.2 Reescribir: el texto describe algo que ya no existe o es falso

- **`ciclo-scf.tex:164-173` (`sec:ci_method`)** describe el CI de parejas nl². Tiene que
  describir la base general de CSFs de pareja: parejas no equivalentes, el pool del Fock de core
  congelado V^(N−2) por cada l ≤ 3 con m = 20, y la autoverificación T1–T4.
  - Además, hay que unificar cómo se llama a la correlación: «estática y dinámica» en `:162`,
    «estática» en `introduccion.tex:22` y «dinámica» en `correlacion-estructura-fina.tex:124`.
- **`ciclo-scf.tex:194-195`:** el propio `\todo` dice que el párrafo no describe lo que se
  hizo. α_d se calibró barriendo contra el NIST, no se tomó de mediciones.
- **`correlacion-estructura-fina.tex:122-140`:** es una exposición genérica de CI completo.
  Tiene que terminar en el modelo que sí se usa y en su límite. Registro a corregir:
  «formidable», «monstruo», «descomunal».
- **`correlacion-estructura-fina.tex:154-167`:** mecanismo de V_pol (B1).
- **`resultados.tex`:** todo lo de la sección 2, más lo siguiente.
  - `:11`: agregar Sn.
  - `:214`: ajustar a la nueva `sec:ci_method`.
  - `:259-311`: toda la discusión del espectro descansa en los números viejos. El ³P₂ de Si pasa
    de +47 % (327.1) a −5.7 % (210.6); el ¹S₀ de C, de 29429.5 a 24574.8; el ¹D₂ de Ge, de
    9648.7 a 7851.6.
  - `:313-317`: rehacer al final.
- **`discusion.tex`** (12-ago) necesita una reescritura completa:
  - `:14` atribuye la desviación de Si al promedio de configuración, pero fue el signo de f₂
    (`HALLAZGOS` §1-2). El origen de esa narrativa está en `docs/avances/estado_proyecto.org:64`.
  - `:19` llama «relatividad escalar» al espín-órbita, que no es un efecto escalar.
  - `:61-69`: el g_eff(Ge) = 1.4975 es de junio y no tiene procedencia. `lande_g_factors.jl:16-17`
    lleva los coeficientes de mezcla a mano, y `docs/figures/lande_g.txt` dice 1.410 para Ge.
    - El log del 07-09 da g_eff(³P₂) = 1.499996 (C), 1.499903 (Si), 1.496164 (Ge) y 1.472668 (Sn).
    - Para Ge, ese valor cae sobre el 1.496 que el texto atribuye al NIST, así que la «brecha»
      de `:69` desaparece.
    - Hay que regenerarlo con procedencia y citar el valor experimental.
  - `:71-73` describe el CI np² → nd².
  - Falta: el cruce de las razones de ζ; que la capa f mueve al ¹D y no al ³P (regla del
    triángulo); y que los virtuales compactos no aportan.
  - Registro: es el archivo con más vocabulario inflado (3.3).
- **`introduccion.tex`** (escribirla al final):
  - `:8` dice «secuencia isoelectrónica», pero C, Si, Ge y Sn son una secuencia homóloga;
    `resultados.tex:4` lo dice bien.
  - `:31-35`: los números de capítulo están escritos a mano y mal. Habla de seis capítulos y
    `main.tex` tiene ocho, más el prólogo sin numerar; conviene usar `\ref{ch:...}`.
  - `:22-27` no menciona el CI completo, C-DIIS ni el resultado de ζ. Además, «un modelo
    robusto que reproduce con alta fidelidad» promete más de lo que da la tabla de ζ.
- **Apéndice F (`metodos-auxiliares.tex`):**
  - `:99-104`: la derivación da W²(r)/r⁴ y dice que es «exactamente» el potencial de
    `eq:vpol_efectivo`, que lleva W. El código usa W (`src/integrals.jl:128-129`,
    `src/atomic_calculations.jl:125-126`). Hay que corregir el texto (por ejemplo, corte √W₆
    sobre el campo). De paso, la ecuación está definida tres veces (`eq:vpol_cutoff`,
    `eq:vpol_efectivo`, `eq:cpp_un_cuerpo`), contra la regla de definir cada ecuación una sola vez.
  - `:119-121` afirma que se usa ARPACK y que BandedMatrices se descartó. No hay `eigs` ni
    `using Arpack` en `src/`, `benchmarks/` ni `examples/`, y `src/AtomicSplines.jl:4` sigue
    haciendo `using BandedMatrices`. Verificar o quitar.

### 3.3 Pasada final

- **`\todo` activos:** 22, más uno comentado en `resultados.tex:98`. Por archivo:
  `galerkin-poisson.tex` 9, `ciclo-scf.tex` 6, `bsplines-galerkin.tex` 3, apéndice `bsplines.tex`
  2, `correlacion-estructura-fina.tex` 1 y `unidades-atomicas.tex` 1.
  - La pregunta de `galerkin-poisson.tex:106` tiene respuesta corta: para una suma finita,
    intercambiar suma e integral es linealidad de la integral; no hace falta convergencia
    uniforme.
- **Citas pendientes:**
  - STO/GTO y condición de Kato (`bsplines-galerkin.tex:4`);
  - la precisión de máquina de Bachau (`:7`);
  - Bachau y el fenómeno de Runge (`:45`);
  - *lifting* (`ciclo-scf.tex:155`);
  - NIST ASD (`ciclo-scf.tex:200`; confirmar la «versión 6.1» de `:205` y poner fecha de consulta);
  - la referencia F.F. del F²(4p,4p) de Ge;
  - el g experimental del ³P₂ de Ge;
  - el origen de α_d = 0.5.
- **Comillas rectas** (`es-noquoting`): `ciclo-scf.tex:58` y `:205`, `galerkin-poisson.tex:81` y
  `:223`, `momento-angular.tex:49`, apéndice `bsplines.tex:118`.
- **Vocabulario inflado.** Es un conteo mecánico de términos como masivo, colapso, contundente,
  formidable, elegante o drástico; hay que revisarlo caso por caso.
  - Por archivo: `discusion.tex` 9, apéndice `bsplines.tex` 8, `correlacion-estructura-fina.tex`
    7, `ciclo-scf.tex` 7, `operadores-tensoriales.tex` 7, `metodos-auxiliares.tex` 5,
    `bsplines-galerkin.tex` 4, `resultados.tex` 3, `particulas-independientes.tex` 3,
    `introduccion.tex` 3 y `momento-angular.tex` 1.
- **`todonotes`** se carga dos veces (`config/preamble.tex:29` y `:39`). Para la versión final
  hay que cargarlo con `[disable]`.
- **Detalles de formato:** `fig:scf-cycle` (`ciclo-scf.tex:38`) nunca se referencia, y la tabla
  de `galerkin-poisson.tex:35` se sale de la página.
- **Petición 3 del asesor (junio), «contrastar con otros métodos»** (`estado_proyecto.org:38`):
  sigue a medias.
  - Solo hay Froese Fischer (límite HF) y NIST.
  - No hay ningún cálculo correlacionado o relativista de la literatura (MCHF-Breit-Pauli, MCDHF,
    CI+MBPT) contra el cual comparar los desdoblamientos.
  - La bibliografía tiene 16 entradas y se citan 13.

---

## 4. Decisiones pendientes (del usuario, antes de regenerar)

1. **Qué orbitales alimentan CI + Breit-Pauli: promedio o ³P.**
   - *Promedio (recomendado).*
     - Los barridos de V_pol ya son de promedio, así que no hace falta ningún SCF nuevo salvo
       para un α_d fuera del barrido.
     - ³P, ¹D y ¹S se tratan con la misma base, que es lo que pide el comentario de
       `np2_toy_ci.jl:430`; `tin_rohf_results_av_R30.0.jld2` se generó justo para esto.
     - Las tablas que validan contra Froese Fischer (energía, ⟨r⁻³⟩, F^k) se quedan en ³P,
       porque comparan el límite HF del término fundamental.
     - Las razones de ζ de trabajo futuro pasan a 1.435, 0.938, 0.884 y 0.857; el cruce se
       conserva.
   - *³P.*
     - HF+CI ya corrió el 07-09.
     - Faltan SCF ³P con V_pol para cada α_d que se reporte: Ge 0.5, y Sn 2–6 si se conserva el
       barrido (cinco SCF de ~6 min con C-DIIS).
     - `tab:estanio_ci` pasa a ³P.
2. **α_d del estaño para la columna CI+V_pol.** `resultados.tex:292` adopta 3.5 «de compromiso»,
   pero no hay SCF con 3.5 y la columna no existe. Conviene decidirlo después de regenerar
   `tab:estanio_ci` con el CI completo: en estaño la mezcla es grande (g_eff = 1.4727) y el
   barrido puede moverse.
3. **Germanio: ¿barrido o solo α_d = 0.5?** También hay que decir de dónde sale el 0.5. Un
   barrido de cinco puntos cuesta del orden de 45 min (SCF de ~2 min y CI de ~6.5 min por punto).
4. **Energía ROHF: autovalor o cociente de Rayleigh.**
   - La cifra de 1.29×10⁻⁵ Ha de `HALLAZGOS` §6 se midió con el signo de f₂ todavía invertido:
     −288.78983… es la energía contaminada, no la canónica −288.85435715. Por eso no sirve para
     decidir.
   - Si la magnitud importa, hay que re-medirla sobre el archivo canónico, lo que tarda segundos.
     Si no importa, se conserva el autovalor y se documenta en el Apéndice F.
   - Va antes de B3: cambiarlo después obliga a correr otra vez los cuatro SCF y todas las tablas.
5. **Formato de la tabla de niveles:** una tabla de cuatro elementos (lo que emite
   `generate_table.jl`, con `\resizebox`) o las dos actuales.
6. **Trabajo futuro:** sección dentro de Conclusiones (recomendado) o capítulo propio.

---

## 5. Orden sugerido

**Fase A — sin cómputo**
1. Hacer `git push` y commitear `CIERRE-DEL-CODIGO.md` y este documento.
2. Tomar las decisiones de la sección 4.
3. Biblioteca:
   - el F²(4p,4p) de Ge en Froese Fischer (1977);
   - el g experimental del ³P₂ de Ge;
   - las citas pendientes;
   - algún cálculo de referencia para la petición 3 del asesor.

**Fase B — lo último que se programa**

4. Escribir una pasada de generación que emita las ocho tablas y las siete figuras desde
   archivos registrados, con la procedencia de cada número. Tiene que cubrir:
   - energía y virial;
   - ⟨r⁻³⟩, F⁰ y F²;
   - Koopmans con y sin CI;
   - niveles HF+CI y +V_pol;
   - `tab:estanio_ci` y la ionización de Ge;
   - g_eff y los mínimos de V_rad;
   - los CSV de convergencia del CI y de C-DIIS para las dos figuras nuevas.
5. Correr B1 sobre los archivos con V_pol ya registrados, con SCF nuevos solo para los α_d que se
   decidan.
6. Correr la pasada de principio a fin, avisando antes, porque tarda del orden de 2–3 h
   dominadas por el estaño.
   - Guardar el log donde git lo vea (hoy `*.log` está ignorado).
   - Pegar el bloque de procedencia junto a cada tabla.
7. Marcar como legado, o borrar, lo que ya no alimenta la tesis:
   - `*_fine_structure.jl`. Su F² a mano y la contradicción de `silicon_fine_structure.jl:163` y
     `:187-195` dejan de importar si nadie los consume.
   - `lande_g_factors.jl`.
   - los `calibrate_*.jl` en su forma actual.
   - `scaling_law.jl`, `plot_energy_levels.jl` y `examples/scratch/plot_*.jl`, una vez que el
     paso 4 los sustituya.

   Opcional: regenerar C y Si `_3P_` para completar `active_s` y `active_p` (nadie las lee).

→ **Congelar el código.**

**Fase C — redacción, en orden de dependencia.** El texto nuevo puede empezar antes de la fase B
usando los números del log del 07-09 como vista previa, siempre que no se peguen en las tablas.

8. `resultados.tex`:
   - tablas y figuras nuevas;
   - los pasajes `:97-98`, `:172-182`, `:206-208` y `:248-317`;
   - subsecciones nuevas de convergencia del CI y de C-DIIS.
9. `ciclo-scf.tex`: `sec:ci_method`, la formulación de C-DIIS, el procedimiento real de
   calibración de V_pol y la cita del NIST.
10. `correlacion-estructura-fina.tex`: el mecanismo de V_pol, la sección de CI y el `\todo`
    de `:4`.
11. `metodos-auxiliares.tex`: la implementación de C-DIIS, W frente a W², y lo de ARPACK y
    BandedMatrices.
12. `discusion.tex`, completo.
13. `conclusiones.tex`, con trabajo futuro.
14. `introduccion.tex` y el resumen.
15. La pasada final de la sección 3.3.

---

## 6. Lo que está bien y no hay que tocar

- **Referencias cruzadas:** ningún `\ref` sin `\label` ni `\label` duplicado; las 13 claves de
  cita existen en `bibliografia.bib` y los 11 archivos de figura existen.
- **Estructura:** cumple las peticiones 1 y 4 del asesor (formato por capítulos; álgebra
  tensorial en apéndices).
- **Datos y motor:** los 16 resultados coinciden con el manifiesto, y el motor de CI se
  autoverifica (T1–T4 a 10⁻¹⁶).
- **Tablas que reproducen:** `tab:momentos_inversos` y la columna ζ de `tab:estanio_ci` salen
  igual desde los archivos registrados.
