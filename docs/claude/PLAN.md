# Plan de trabajo: CI de valencia completo + C-DIIS

Escrito 2026-09-06 tras una auditoría del motor de CI. Léelo entero antes de tocar nada.
Las reglas duras y el mapa del repo están en `CLAUDE.md` (raíz).

---

## 1. De dónde venimos

La tesis calcula el espectro de multipletes de la secuencia $np^2$ (C, Si, Ge, Sn): ROHF con
B-splines → integrales de Slater $F^k$ y parámetro espín-órbita $\zeta_{np}$ → CI de valencia →
diagonalización de Breit-Pauli en acoplamiento intermedio.

Una auditoría encontró que **el CI no estaba haciendo nada**. Recuperaba entre 170 y 1036 cm⁻¹ de
correlación (C: $-7.75\times10^{-4}$ Ha, Si: $-4.72\times10^{-3}$, Ge: $-2.76\times10^{-3}$), no
monótono en $Z$, y uno o dos órdenes por debajo de la correlación de valencia real (~0.04–0.07 Ha).
El factor de Landé efectivo con CI y sin CI coincidía en la séptima cifra. **La columna "HF+CI" de
`tab:niveles_energia_ligeros` es, número por número, la columna HF.**

La causa era el espacio de configuraciones: solo admitía CSFs de pareja *equivalente* $(nl)^2$, con
2 orbitales $s$, 3 $p$ y 2 $d$, y sin CSFs no equivalentes ($ns\,np$, $np\,n'p$, $nd\,n'd$) ni $l=3$.
Con una base de splines, que da un pseudo-espectro completo por cada $l$, usar dos o tres orbitales
es tirar justo lo que la base regala.

## 2. Qué ya está hecho y verificado

`benchmarks/np2_sequence/np2_toy_ci.jl` (respaldo en `.bak`) fue parcheado y **el usuario ya lo corrió**:

- Fase angular corregida a $(-1)^{l+l'+L}$ (antes $(-1)^L$; inocuo para autovalores, pero incorrecto).
- $\zeta$ se calcula del `.jld2` en lugar de venir como literal. Método nuevo de dos argumentos
  `run_toy_ci(element, file)`; el de tres se conservó para no romper `calibrate_vpol.jl` y
  `calibrate_tin_vpol.jl`, que sí lo calculaban bien.
- Espacio activo reortonormalizado en la métrica $S$ (el orbital SCF de valencia se insertaba entre
  autovectores del Fock de core congelado y rompía la ortonormalidad que el CI asume).
- La raíz física se elige por máximo $|C_{np^2}|$, no por ser la más baja del bloque.
- **Test de consistencia angular**: sobre el CSF de referencia, $[E(^1D)-E(^3P)]/F^2$ debe valer
  exactamente 0.24 y $[E(^1S)-E(^3P)]/F^2$ exactamente 0.60.

**Resultado de esa corrida (estaño):** deltas de $-3.9\times10^{-16}$ y $-2.2\times10^{-16}$, y
$F^2(5p,5p) = 0.14723221$ Ha frente al $0.14723200$ del Cuadro `tab:integrales_slater`.
**`compute_Rk` y los coeficientes angulares están sanos.** Una sospecha previa de un déficit del 4–6%
en $F^2$ resultó ser un artefacto de un log viejo; queda descartada.

`generate_table.jl` también se reescribió: llamaba a una firma inexistente y fallaba con `MethodError`,
es decir, **las tablas del capítulo 6 no se podían regenerar**. Ahora lleva la procedencia de cada
archivo escrita en el propio script y aborta si el test angular falla.

## 3. Qué está escrito pero NUNCA se ha ejecutado

`benchmarks/np2_sequence/np2_ci_full.jl` (~450 líneas). CI de pareja de valencia general en $(n,l)$:

- `PairCSF(i,j,L,S)` con $i \le j$ sobre un pool de orbitales; `generate_pair_csfs` aplica paridad par,
  regla del triángulo y Pauli ($L+S$ par para electrones equivalentes).
- Normas $N_{ab}=1/2$ (equivalentes) o $1/\sqrt2$; fase de intercambio $(-1)^{l_a+l_b+L+S}$.
- Elemento de matriz
  $2N_{ab}N_{cd}\left[\langle ab|H|cd\rangle + (-1)^{l_c+l_d+L+S}\langle ab|H|dc\rangle\right]$,
  con el término directo
  $\sum_k (-1)^{l_b+l_c+L}\begin{Bmatrix} l_a & l_b & L \\ l_d & l_c & k\end{Bmatrix}
   \langle l_a\|C^k\|l_c\rangle\langle l_b\|C^k\|l_d\rangle R^k(ab,cd)$.
- `build_orbital_pool` diagonaliza el Fock de core congelado $V^{N-2}$ por cada $l\le3$, empalma el
  orbital SCF de valencia y reortonormaliza.
- `engine_selftest` con tres pruebas, y aborta si alguna falla:
  - **T1** hermiticidad de cada bloque.
  - **T2** sobre el subespacio de parejas equivalentes el motor general debe reproducir *dígito a
    dígito* al motor legado (`interaction_coefficient`). Valida fases, normas y 6j contra código ya
    verificado. **Es la prueba que importa.**
  - **T3** los coeficientes exactos de $p^2$ (0.24 / 0.60 $F^2$).
- `convergence_study` sube el espacio activo y tabula $E_{corr}$ frente al tamaño.

**Límite de costo, conocido y aceptado:** `compute_Rk` resuelve una ecuación de Poisson por integral,
sin lista transformada de integrales bi-electrónicas. El número de $R^k$ distintos crece como
$N_{orb}^4/\text{simetría}$, así que el pseudo-espectro completo **no es alcanzable** sin reescribir
el motor de integrales. El entregable no es "el CI da X" sino "$E_{corr}$ converge a X con este
espacio activo, y aquí está la curva". Si hiciera falta ir más lejos, el siguiente paso de ingeniería
es precomputar y almacenar la lista de $R^k$ una sola vez, no optimizar el bucle.

---

## 4. Fases

### Fase 0 — Compuerta: que el motor nuevo se valide a sí mismo

```bash
julia --project=. benchmarks/np2_sequence/np2_ci_full.jl 2>&1 | tee /tmp/ci_full_carbon.log
```

Arranca con carbono (100 elementos, $K=7$: el más barato), `sizes=[3,4,6]`, `lmax=2`.

**No sigas hasta que T2 pase.** Si falla, el error casi seguro está en una de estas tres cosas, en
este orden de probabilidad: el orden de argumentos del 6j en `g_direct`
(`wigner6j(la, lb, L, ld, lc, k)` debe ser $\{l_a\,l_b\,L;\,l_d\,l_c\,k\}$); la fase
$(-1)^{l_b+l_c+L}$ combinada con el $(-1)^{l_a+l_b}$ de los elementos de matriz reducidos; o el
prefactor $2N_{ab}N_{cd}$. Depura comparando un solo elemento equivalente a mano.

Si T1 falla pero T2 pasa, el problema es de simetría en `h_direct`, no en la parte angular.

### Fase 1 — Convergencia del CI

Con T2 verde, corre el estudio de convergencia y responde: **¿$E_{corr}$ sube a valores del orden de
$10^{-2}$ Ha y se estabiliza?**

- Sube `lmax` a 3 y el tamaño hasta donde aguante el tiempo de cómputo. Anota la curva completa,
  incluyendo el número de $R^k$ y el tiempo de pared: eso va a la tesis como evidencia de convergencia.
- Vigila `pick_reference_root`: si empieza a avisar que la raíz de mayor peso no es la más baja,
  es información física real (la base activa está metiendo estados $s^2$ o Rydberg por debajo),
  no un bug. Anótalo.
- Repite para Si, Ge y Sn. El estaño es el caro (500 elementos, $K=8$).

**Preguntas abiertas que esta fase debería cerrar:**

- **¿Cuál es el $\zeta$ correcto?** `carbon_fine_structure.jl` usa $\zeta(2p) = 31.946$ cm⁻¹
  (Froese Fischer); el CI calculaba $0.000189$ Ha $= 41.47$ cm⁻¹ con `compute_zeta`. Con 31.946 el
  $^3P_2$ del carbono daría ~47.9 cm⁻¹ en vez de 62.2 (NIST: 43.4). Hay que decidir cuál es el bueno
  y por qué, no promediarlos.
- **Auto-interacción en $\zeta$.** `compute_effective_central_potential` se llama con la lista
  **completa** de orbitales, valencia incluida (`carbon_rohf.jl:181`, `tin_rohf.jl:248`, etc.). El
  $V_{eff}$ que alimenta el operador espín-órbita incluye el apantallamiento del electrón de valencia
  sobre sí mismo; debería ser el potencial de core ($N-1$). **Esto se puede arreglar sin volver a
  correr el SCF**: los `orbitals` están en el `.jld2`, basta recomputar `V_eff` quitando la valencia.
  Corregirlo *sube* $\zeta$, así que no explica la sobrestimación, pero es un error metodológico real
  y hay que cuantificarlo.

### Fase 2 — C-DIIS en el ROHF (independiente de las fases 0 y 1, puede ir en paralelo)

**Motivación honesta para la tesis:** el ROHF de Ge y Sn tarda ~200 iteraciones. La razón está a la
vista en `tin_rohf.jl:118`: `level_shift = 3.0` se aplica en **cada** iteración y **nunca se apaga**.
Un desplazamiento constante de 3 Ha sobre el espacio virtual amortigua la convergencia de forma
permanente: es una muleta que nunca se suelta. La corrección es soltar la muleta cuando el sistema ya
está en la cuenca de atracción y cambiar a extrapolación de Pulay.

**Prototipo existente que funciona:** `benchmarks/closed_shell/11_radon.jl`, líneas ~112–215.
Vector de error conmutador $\mathbf{e} = FDS - SDF$ por bloque de $l$, matriz de Pulay $B$ con
productos de Frobenius sumados sobre bloques, sistema aumentado con multiplicador de Lagrange,
extrapolación $F_{\text{eff}} = \sum_i c_i F_i$, historial de 6, y sin mezcla una vez que DIIS entra.
Ese código es de **capa cerrada**: hay que adaptarlo, no copiarlo.

**Qué cambia en el caso ROHF de la secuencia $np^2$:**

`tin_rohf.jl` construye **cuatro** operadores de Fock, no uno por $l$:

| operador | actúa sobre | construido con |
|---|---|---|
| `F_s` | 1s…5s (cerradas) | `J_core + J_5p_spherical - K_core_s - K_5p_on_s` |
| `F_core_p` | 2p, 3p, 4p (cerradas) | `J_core + J_5p_spherical - K_core_p - K_5p_on_p` |
| `F_core_d` | 3d, 4d (cerradas) | `J_core + J_5p_spherical - K_core_d - K_5p_on_d` |
| `F_5p` | 5p (**abierta**, occ 2 de 6) | `J_core + J_5p_intra - K_core_p - K_5p_intra` |

De ahí las reglas de la adaptación:

1. **Cada matriz de error se construye con su pareja $(F, D)$ correcta.** $D_{5p}$ sale del orbital de
   valencia con su ocupación fraccionaria; $D_{core\_p}$ de 2p/3p/4p. Mezclarlas invalida el
   conmutador.
2. **El error y la extrapolación se hacen sobre el Fock SIN desplazar.** El level shift es una muleta
   numérica, no parte del operador convergido, y sus proyectores cambian en cada iteración: extrapolar
   Focks desplazados de iteraciones distintas es inconsistente. Guarda el $F$ crudo en el historial y
   aplica el shift, si toca, después de extrapolar.
3. **Criterio de cambio de régimen.** Mantén level shift y SCF simple hasta que
   $\max_\text{bloques}\|\mathbf{e}\|_\infty$ baje de un umbral (empieza con $10^{-2}$), y entonces
   apaga el shift y enciende DIIS. No lo hagas por número fijo de iteraciones: hazlo por residual, y
   registra en qué iteración ocurrió.
4. **Salvaguardas.** Si $\text{cond}(B)$ se dispara, tira los vectores más viejos; si el solve falla o
   $\sum_i c_i$ se aleja de 1, cae a una iteración de SCF simple y reinicia el historial. Sin esto,
   DIIS diverge silenciosamente y es peor que no tenerlo.
5. **Advertencia física que hay que decir en la tesis:** que $FDS-SDF$ se anule bloque por bloque es
   condición necesaria de estacionariedad dentro de cada canal $l$, pero **no captura el acoplamiento
   entre la capa abierta $5p$ y las cerradas**. Esto es "C-DIIS por bloques sobre operadores de Fock
   dependientes del término", no C-DIIS canónico de Roothaan. Funciona en la práctica, pero la
   convergencia se declara con la energía y el cociente virial, no con el residual de DIIS.

**Entregable de la fase:** una tabla de iteraciones hasta $|\Delta E| < 10^{-10}$ Ha para C, Si, Ge y
Sn, con y sin C-DIIS, a la misma tolerancia y desde el mismo punto de partida. **Comprueba que ambos
métodos llegan a la misma energía en las 8 cifras que la tesis ya reporta contra Froese Fischer**: si
C-DIIS converge más rápido a un número distinto, converge a otra solución y eso es un hallazgo, no una
mejora.

### Fase 3 — Recalibración de $V_{pol}$

Con el CI bueno y el SCF rápido, rehacer el barrido de $\alpha_d$ para Ge y Sn
(`calibrate_vpol.jl`, `calibrate_tin_vpol.jl`), que ahora deben llamar al CI completo en vez del de
pares. Cada valor de $\alpha_d$ necesita su propio SCF: **es aquí donde C-DIIS se paga solo**, seis o
más SCF completos por elemento.

### Fase 4 — Tablas y manuscrito

Regenerar con `generate_table.jl` (rellenando las rutas `file_vpol`, hoy en `nothing`) y pegar el
bloque de procedencia como comentario junto a cada tabla en el `.tex`.

Cambian `tab:niveles_energia_ligeros`, `tab:niveles_energia_pesados`, `tab:estanio_ci` y
posiblemente `tab:integrales_slater` (cuyo renglón de Germanio, además, es un copy-paste dígito a
dígito del de Silicio: 0.16272437 en ambos). Hay que reescribir al menos el párrafo de
`resultados.tex:311` (atribuye el 3.1% residual al espacio CI truncado), la discusión de $V_{pol}$ en
`resultados.tex:263–269`, y `sec:ci_method` de `ciclo-scf.tex`, que describe un CI de pares y debe
pasar a describir la base de CSFs no equivalentes.

---

## 5. ¿Hay que volver a correr los SCF?

- **Para el trabajo de CI: no.** Los `orbitals` convergidos están en los `*_rohf_results_*.jld2`.
- **Para arreglar la auto-interacción en $\zeta$: no.** `V_eff` se recomputa post hoc desde los
  `orbitals` guardados, quitando la valencia.
- **Para la recalibración de $V_{pol}$ (fase 3): sí**, un SCF por cada valor de $\alpha_d$.
- **Para la tabla de convergencia de C-DIIS (fase 2): sí**, pero ese *es* el experimento.

O sea: el motivo real para meter C-DIIS no es desatascar el trabajo de CI, es que la fase 3 son
docenas de SCF y hoy cada uno cuesta ~200 iteraciones.

## 6. Dónde va C-DIIS en la tesis

Siguiendo la regla editorial "el capítulo muestra los resultados, el apéndice los pasos":

- **Formulación** (vector de error de Pulay, sistema aumentado, por qué el conmutador mide
  estacionariedad) en `docs/chapters/ciclo-scf.tex`. Es parte legítima del ciclo autoconsistente y el
  capítulo tiene espacio; ahí ya vive la estructura algebraica del SCF.
- **Resultado numérico** (tabla de iteraciones con y sin C-DIIS, y la curva de residual) en
  `resultados.tex`, con la advertencia del punto 5 de la fase 2.
- **Detalles de implementación** (tamaño del historial, condicionamiento de $B$, criterio de cambio de
  régimen, salvaguardas) en el Apéndice F, `metodos-auxiliares.tex`.

No inventes la sección: escríbela solo cuando existan los números de la fase 2.
