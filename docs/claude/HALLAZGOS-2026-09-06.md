# Hallazgos del 2026-09-06

Sesión de trabajo sobre las fases 0, 2 y 1 de `PLAN.md`, en ese orden. Este documento
recoge lo que se encontró y lo que queda por decidir; no reemplaza al plan.

---

## 1. El hallazgo que más pesa: el signo de $f_2$ en carbono y silicio

`carbon_rohf.jl` y `silicon_rohf.jl` tenían `coeff_k2 = -5/25` para el estado $^3P$,
mientras que `germanium_rohf.jl` y `tin_rohf.jl` usaban `+5/25`, con el mismo ensamblado
(`F_np .= ... .- K_np_intra`) y la misma llamada. **Germanio y estaño eran los correctos.**

El álgebra fija la convención sin ambigüedad. Slater-Condon da $f_2(^3P) = -5/25$, que es
justo lo que verifica el test T3 del motor de CI (las diferencias $^1D-^3P$ y $^1S-^3P$
valen 0.24 y 0.60 $F^2$). Como el ensamblado resta $K_{intra}$, la convención es
$\text{coeff\_k2} = -f_2$, y el promedio de configuración pesado por degeneración,
$\langle f_2\rangle = [9(-5)+5(1)+1(10)]/375 = -2/25$, es coherente con el `+2/25` que los
cuatro scripts ya usaban para `av`. Luego el $^3P$ exige `+5/25`.

Tres pruebas independientes lo confirman:

1. **Ordenamiento imposible.** $E(^3P) - E(\text{av})$ salía positiva en carbono
   ($+0.0653$ Ha) y silicio ($+0.0448$ Ha): el término fundamental por encima del
   promedio de configuración. Ambas caen sobre la predicción del signo invertido,
   $+0.28\,F^2$, con un 2% de diferencia atribuible a relajación orbital. Germanio da
   $-0.0194$ contra la predicción correcta $-0.12\,F^2 = -0.0195$.

2. **Límite Hartree-Fock.** Con el signo corregido, carbono da $E = -37.68861896$ Ha y
   silicio $-288.85435715$ Ha, que son los valores de Froese Fischer **y** los que ya
   reporta `resultados.tex:25` y `:30`. Con el signo viejo el código daba $-37.594380$ y
   $-288.789838$, valores que no aparecen en ninguna parte del manuscrito. El error
   variacional era de 0.094 Ha en carbono y 0.065 Ha en silicio.

3. **Integrales de Slater.** El $F^2(3p,3p)$ corregido de silicio vale $0.16594819$, a
   $1.6\times10^{-5}$ de Froese Fischer; el valor que hoy aparece en la tesis
   ($0.16272437$) está a $3.2\times10^{-3}$.

4. **Corroboración independiente.** `carbon_fine_structure.jl:128` lleva escrito a mano
   `F2_val = 0.24330175`, que es exactamente el $F^2$ que produce el código corregido y el
   que aparece en `tab:integrales_slater`. Ese script se escribió cuando el código estaba
   bien. (De paso: es una cantidad física congelada que el código puede calcular, o sea un
   caso de la regla 3 de `CLAUDE.md`, aunque el valor sea correcto.)

**Consecuencia:** el código regresionó *después* de generar las tablas del capítulo 6. Los
archivos `carbon_rohf_results_3P_R30.0.jld2` y `silicon_rohf_results_3P_R30.0.jld2` se
regeneraron con el código roto y no son consistentes con el manuscrito. Todo lo que salga
de ellos —$\zeta$, $F^2$, el CI completo— hereda el error. Germanio y estaño no están
afectados, y los archivos `_av_` de los cuatro tampoco, porque `+2/25` siempre fue correcto.

### Estado

- Corregido en `carbon_rohf.jl:40` y `silicon_rohf.jl:39`, con el razonamiento escrito al lado.
- SCF regenerados en archivos **nuevos**: `carbon_rohf_results_3P_R30.0_f2fix.jld2` y
  `silicon_rohf_results_3P_R30.0_f2fix.jld2`. Los originales no se tocaron, porque los
  `.jld2` están en `.gitignore` y sobrescribirlos sería irreversible.
- Falta decidir si se reemplazan los originales.

---

## 2. Procedencia rota en `tab:integrales_slater`

| | en la tesis | recalculado | Froese Fischer |
|---|---|---|---|
| Si $F^0(3p,3p)$ | 0.31650266 | **0.32969971** | 0.32971786 |
| Si $F^2(3p,3p)$ | 0.16272437 | **0.16594819** | 0.16593265 |
| Ge $F^0(4p,4p)$ | 0.30137782 | **0.31650266** | 0.32624597 |
| Ge $F^2(4p,4p)$ | 0.16272437 | 0.16272437 | 0.16593265 |

Los valores del renglón de **silicio** son, dígito a dígito, los que el código calcula
para **germanio**. El copy-paste va de Ge a Si, al revés de lo que suponía `PLAN.md`. El
$F^0$ del renglón de germanio (0.30137782) no reproduce con ningún cálculo: procedencia
desconocida, hay que rastrearlo antes de volver a publicarlo. Carbono y estaño sí
reproducen, a $10^{-8}$ y $10^{-5}$ respectivamente.

La discrepancia de silicio que el texto atribuye a la discretización es en realidad el
signo de $f_2$ más una transcripción equivocada. Corregidas ambas, silicio cae a
$1.6\times10^{-5}$ de la referencia.

**La columna de Froese Fischer también está duplicada.** El $F^2$ de referencia de silicio
y el de germanio son el mismo número a ocho cifras (0.16593265), lo cual no puede ser para
orbitales $3p$ y $4p$ de átomos distintos; el $F^0$ de referencia sí difiere entre ambos.
Hay que volver a la fuente y recuperar el $F^2(4p,4p)$ de germanio antes de republicar la
tabla. Con la referencia duplicada, el $\Delta$ de $3.2\times10^{-3}$ que hoy aparece en
el renglón de germanio no significa nada.

Los cuatro valores calculados con orbitales correctos, para cuando se rehaga la tabla:

| | $F^0$ | $F^2$ |
|---|---|---|
| C  | 0.53860369 | 0.24330175 |
| Si | 0.32969971 | 0.16594819 |
| Ge | 0.31650266 | 0.16272437 |
| Sn | 0.27874830 | 0.14723221 |

---

## 3. Fase 0: el motor de CI está verificado, la autoverificación no lo estaba

`np2_ci_full.jl` corrió a la primera y **T2 pasó a precisión de máquina** ($5.6\times10^{-17}$).
Las tres sospechas del plan sobre el orden del 6j, la fase y el prefactor quedan descartadas.

Lo que sí estaba mal era la autoverificación:

- **T1 no podía fallar.** `build_pair_hamiltonian` llenaba `H[i,j]` y `H[j,i]` con la misma
  llamada, así que `max|H-H'|` era cero por construcción y la rama "si T1 falla pero T2
  pasa" del plan era inalcanzable. Ahora reevalúa el triángulo superior con los CSFs
  invertidos y mide $4.7\times10^{-16}$.
- **Faltaba validar lo único que el motor agrega.** T2 solo cubre pares equivalentes, que
  es exactamente lo que ya hacía el motor legado. Se añadió **T4**, que contrasta pares no
  equivalentes contra Condon-Shortley: para $(ns, n'l')$ debe cumplirse
  $E(^1L)-E(^3L) = 2G^{l'}/(2l'+1)$. Valida normas $1/\sqrt2$, fase de intercambio y 6j con
  $l_a \neq l_b$. Pasa a $10^{-16}$ en `3s4s`, `3s2p` (paridad impar), `3s3d` y `3s4f`.

---

## 4. Fase 1: convergencia del CI

`convergence_study` ahora comparte el workspace entre tamaños (los pools son anidados y la
clave del caché es $(n,l,k)$, así que el espacio grande no recalcula los $R^k$ del chico) y
reporta orbitales, CSFs, $R^k$, tiempo de pared y las raíces elegidas.

**Respuesta a la pregunta del plan: sí, $E_{corr}$ llega al orden de $10^{-2}$ Ha y converge.**

Resultado completo a $m = 20$, `lmax = 3` (80 orbitales, 630 CSFs en el bloque $^3P$,
$\sim 2.5\times10^6$ integrales $R^k$). C y Si sobre orbitales corregidos:

| | $E_{corr}$ (Ha) | $^3P_1$ | $^3P_2$ | $^1D_2$ | $^1S_0$ | $t$ (s) |
|---|---|---|---|---|---|---|
| C  calculado | $-0.00870$ |   21.15 |   63.08 | 11148.4 | 24574.8 |  92 |
| C  NIST      |            |   16.40 |   43.40 | 10192.6 | 21648.0 |     |
| Si calculado | $-0.00856$ |   72.35 |  210.58 |  7285.2 | 14674.5 |  83 |
| Si NIST      |            |   77.10 |  223.20 |  6298.8 | 15394.4 |     |
| Ge calculado | $-0.00706$ |  490.74 | 1257.58 |  7851.6 | 15780.2 | 390 |
| Ge NIST      |            |  557.10 | 1410.00 |  7125.3 | 16367.1 |     |
| Sn calculado | $-0.00624$ | 1367.10 | 2921.70 |  8710.1 | 16045.3 | 609 |
| Sn NIST      |            | 1691.80 | 3427.70 |  8613.0 | 17162.6 |     |

La calidad de la convergencia **crece con $Z$**. Sumando la cola geométrica de los
incrementos, lo que falta por recuperar es $-1.5\times10^{-3}$ Ha en carbono,
$-8.8\times10^{-5}$ en silicio, $-3.8\times10^{-5}$ en germanio y $-3.1\times10^{-5}$ en
estaño: los dos pesados ya están convergidos a $m=20$, carbono no. En los cuatro casos y
en todos los tamaños `pick_reference_root` devolvió 1/1/1, es decir que la raíz de mayor
peso $np^2$ nunca dejó de ser la más baja del bloque; el aviso que anticipaba el plan no
llegó a dispararse ni con 80 orbitales.

Tres patrones que hay que discutir en el capítulo:

1. **La estructura fina sigue la razón de $\zeta$ y nada más**: carbono $+45\%$, luego
   $-5.7\%$, $-10.8\%$ y $-14.8\%$. Es el mismo cruce monótono de la sección 5.
2. **El $^1D_2$ queda alto en toda la secuencia** ($+9.4\%$, $+15.7\%$, $+10.2\%$), salvo
   estaño, que da $+1.1\%$.
3. **El $^1S_0$ cambia de signo**: alto solo en carbono ($+13.5\%$), bajo en los otros tres.

Dos cosas que hay que decir sobre el espacio activo:

- **La capa $f$ importa para $^1D$ pero no para $^3P$.** Aporta $-6\times10^{-6}$ Ha a
  $E_{corr}(^3P)$ y en cambio mueve $^1D_2$ en $-269$ cm$^{-1}$. La razón es la regla del
  triángulo: los CSFs $(p,f)$ están permitidos en el bloque $^1D$ y prohibidos en el $^3P$.
- **Los virtuales compactos de alta energía no aportan nada.** El acoplamiento va como la
  densidad de transición $P_{np}P_v$, que un estado muy oscilante cancela aunque esté
  localizado donde el orbital de valencia. Todo el peso está en los virtuales más bajos.

### Límite estructural del modelo

El CI correla **solo la pareja de valencia** con el core congelado, así que por construcción
no contiene la casi-degeneración $ns^2 \leftrightarrow np^2$, que es un problema de cuatro
electrones y la contribución más grande a la correlación de valencia de la segunda fila.
Como $ns^2$ es $^1S$, esa omisión deprime específicamente al $^1S$.

Esto importa para `resultados.tex:311`, que atribuye el residuo al espacio CI truncado. Hay
que reescribirlo: parte del residuo es del modelo, no del tamaño del espacio, y levantarlo
exige un CI multi-referencia con el $ns$ activo.

---

## 5. La estructura fina es un problema de $\zeta$, no del CI

$^3P_2$ se mueve de 63.15 a 63.09 cm$^{-1}$ en carbono mientras el espacio activo se
quintuplica. El CI no toca la estructura fina.

Alimentando Breit-Pauli con los desdoblamientos LS del NIST, para aislar el álgebra de los
errores del CI, el $\zeta$ que reproduce el $^3P_2$ medido es:

| | Z | $\zeta$ calculado | $\zeta$ que exige NIST | razón |
|---|---|---|---|---|
| C  |  6 |   42.216 |   28.909 | 1.460 |
| Si | 14 |  140.992 |  148.093 | 0.952 |
| Ge | 32 |  823.671 |  917.140 | 0.898 |
| Sn | 50 | 1885.936 | 2168.759 | 0.870 |

(C y Si con los orbitales corregidos.) **No hay un error sistemático único:** la razón cae
monótonamente con $Z$ y cruza 1 entre carbono y silicio. Eso apunta a dos omisiones
distintas de la truncación de Breit-Pauli, cada una dominante en un extremo. En átomos
ligeros faltan los términos de dos cuerpos (espín-otra-órbita y espín-espín), que *reducen*
el desdoblamiento y pesan como $\sim 1/Z_{ef}$ frente al de un cuerpo. En los pesados falta
la contracción relativista del orbital de valencia cerca del núcleo, que lo *aumenta* y
crece con $Z$.

Descartado como causa, todo medido:

- **Auto-interacción en $V_{eff}$**: $+1.8\%$ en carbono, $<0.2\%$ en el resto. Es un error
  metodológico real (el $V_{eff}$ que alimenta al operador espín-órbita incluye el
  apantallamiento del electrón de valencia sobre sí mismo) pero irrelevante numéricamente,
  y sube $\zeta$, que ya sobra en carbono. El control reproduce el `V_eff` guardado bit a bit.
- **La malla de `compute_zeta`**: converge a 6 cifras, insensible a `rmin` y al número de puntos.
- **La elección de orbital $^3P$ frente a promedio de configuración**: cambia $\zeta$ un ~4%.

**Decisión tomada:** la tesis usa el $\zeta$ calculado y explica la discrepancia como
deficiencia del modelo en una sección de trabajo futuro. No se sustituye por el valor de
Froese Fischer ni por uno ajustado al experimento.

---

## 6. Fase 2: C-DIIS

Módulo `rohf_diis.jl` y banco `diis_benchmark.jl`, con el gancho detrás de `use_diis` en
los cuatro `*_rohf.jl` sin tocar el ensamblado de las Fock. El banco corre con
`save = false` y no pisa los resultados.

| | iter sin | iter con | factor | s sin | s con |
|---|---|---|---|---|---|
| C | 26 | 12 | 2.2× | 8.7 | 0.4 |
| Si | 27 | 28 | 1.0× | 3.5 | 1.7 |
| Ge | 241 | 85 | 2.8× | 340.7 | 128.0 |
| Sn | 250 | 52 | 4.8× | 1593.5 | 347.6 |

Los renglones de C y Si se volvieron a medir **después** de corregir el signo de $f_2$: con
el operador roto salían 45→13 y 40→44, porque el Fock equivocado también estaba peor
condicionado y la línea base necesitaba más iteraciones. Ge y Sn no cambian, y son los que
llevan el argumento.

Los dos brazos coinciden entre sí a $\le 2.7\times10^{-9}$ Ha y con los `.jld2` a
$\le 1.5\times10^{-8}$: 11 cifras significativas, muy por encima de las 8 que exigía el
criterio del plan. Silicio queda en tablas porque no tiene level shift que soltar.

Tres cosas que el plan no anticipaba:

1. **El conmutador crudo $FDS-SDF$ no se anula en el punto fijo de este SCF.** El orbital
   de valencia se reortogonaliza por Gram-Schmidt contra las capas cerradas del mismo $l$,
   así que en la solución cumple $F_{np}|np\rangle = \varepsilon|np\rangle + \sum_c \lambda_c |c\rangle$.
   Los multiplicadores de Lagrange le dejan un piso no nulo, medido en 1.2e-4 (Si),
   2.2e-6 (Ge) y 7.1e-7 (Sn). Carbono es el control: es el único sin Gram-Schmidt en el
   canal de valencia y el único sin piso (2.5e-13). Proyectando el conmutador fuera del
   espacio ocupado del canal, el término de multiplicadores muere y queda el bloque
   ocupado-virtual, que sí es cero exacto en la solución ROHF.
2. **La energía usaba autovalores de la Fock extrapolada**, así que dejaba de ser un
   funcional de los orbitales y heredaba el temblor de la extrapolación. Corregirla al
   cociente de Rayleigh lo arregla, pero mueve la línea base de silicio $1.29\times10^{-5}$ Ha
   respecto al manuscrito, así que se revirtió. En su lugar, el test de convergencia
   rechaza $|\Delta E| < \text{tol}$ mientras el residual esté estancado: sin eso silicio
   "convergía" en falso sobre un ciclo límite de periodo 4.
3. **DIIS se congela** cuando su error toca el piso (iteración 21 en Ge y Sn) y el tramo
   final lo hace el SCF simple. La ganancia se reparte entre acelerar el arranque y soltar
   el level shift para la cola.

Las trazas de convergencia están en `diis_trace_<elemento>_{con,sin}_diis.csv`.

---

## 7. Qué queda abierto

1. **Reemplazar o no** los `.jld2` de C y Si por los `_f2fix`. Los originales no se tocaron.
2. **Rastrear el $F^0$ de germanio** de `tab:integrales_slater` (0.30137782), que no
   reproduce con ningún cálculo.
3. **`np2_toy_ci.jl` ejecuta los archivos `_3P_`** pero su propio comentario dice que el CI
   se corre mejor sobre los de promedio de configuración. El comentario y el código se
   contradicen, y no existe `tin_rohf_results_av_R30.0.jld2`, que es por qué la columna del
   estaño de `tab:niveles_energia_pesados` está vacía.
4. **La fila "HF + CI (Valencia)" de la tabla de ionización de germanio** (`resultados.tex:305`)
   da 7.40 eV contra 7.33 de HF puro, o sea 0.0026 Ha de correlación: exactamente lo que
   recuperaba el CI roto. Esa tabla y su párrafo van a cambiar.
5. **`generate_table.jl`** ya llama al CI completo con `ACTIVE_SPACE` documentado, pero
   regenerar la tabla entera cuesta del orden de 40 minutos, dominados por el estaño.
