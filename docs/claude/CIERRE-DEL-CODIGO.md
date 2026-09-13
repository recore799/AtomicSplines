# Cierre del código: qué falta para poder escribir

Escrito 2026-09-08, después de las fases 0, 1 y 2. Este documento define **cuándo se congela el
código**. Todo lo que no esté en la lista de bloqueantes es redacción, no programación.

## Lo que ya está resuelto y no se vuelve a tocar

- **Motor de CI verificado.** T1 (hermiticidad), T2 (reproduce al motor legado sobre parejas
  equivalentes), T3 (0.24/0.60 $F^2$) y T4 (Condon-Shortley para parejas no equivalentes) pasan en los
  cuatro elementos, con residuos de $10^{-16}$.
- **CI convergido.** A $m=20$, `lmax=3` (80 orbitales, 630 CSFs en $^3P$, $2.5\times10^6$ integrales
  $R^k$): $E_{corr} = -8.70, -8.56, -7.06, -6.24 \times 10^{-3}$ Ha para C, Si, Ge, Sn. La cola
  geométrica que falta es $-8.8\times10^{-5}$ (Si), $-3.8\times10^{-5}$ (Ge), $-3.1\times10^{-5}$ (Sn):
  **convergido**. Carbono es el único con cola apreciable, $-1.5\times10^{-3}$ Ha.
- **Signo de $f_2$** en `carbon_rohf.jl` y `silicon_rohf.jl` corregido; el HF vuelve a reproducir a
  Froese Fischer y a lo que ya dice `resultados.tex:25` y `:30`.
- **C-DIIS** funcionando: 2.8× en Ge, 4.8× en Sn, con los dos brazos coincidiendo a $10^{-9}$ Ha.

## Decisiones de alcance tomadas (NO reabrir)

1. **El CI se queda en correlación de valencia con core congelado.** Está convergido respecto al
   espacio activo; lo que falta no es espacio sino modelo. Ampliar el espacio activo no mueve la
   estructura fina ($^3P_2$ pasa de 63.15 a 63.09 cm⁻¹ mientras el espacio se quintuplica).
2. **No se implementa correlación de core.** Justificación física, no de conveniencia: las cantidades
   que reporta la tesis son diferencias de energía dentro de una misma configuración, y la correlación
   core-core es prácticamente idéntica para $^3P$, $^1D$ y $^1S$, de modo que se cancela en los
   desdoblamientos. La parte que **no** cancela es la core-valencia, y de eso ya se hace cargo
   $V_{pol}$ como modelo fenomenológico.
3. **Se usa el $\zeta$ calculado**, no el de Froese Fischer ni uno ajustado. La discrepancia se explica
   como truncación de Breit-Pauli.
4. **La casi-degeneración $ns^2 \leftrightarrow np^2$ queda fuera** y se declara explícitamente como
   límite del modelo: es un problema multi-referencia de cuatro electrones. Deprime específicamente al
   $^1S$, lo cual es consistente con el patrón observado.

## Bloqueantes (esto sí hay que hacerlo antes de escribir el capítulo 6)

### B1. La narrativa de $V_{pol}$ está invertida — PRIORIDAD MÁXIMA

`resultados.tex:263-265` y el caption de `tab:niveles_energia_pesados` dicen que HF **sobrestima** el
desdoblamiento del germanio por un factor ~1.5 y que $V_{pol}$ lo **baja** relajando el orbital hacia
afuera. Los números nuevos dicen lo contrario:

| Ge | $^3P_1$ | $^3P_2$ |
|---|---|---|
| HF+CI (ζ calculado) | 490.7 | 1257.6 |
| CI+$V_{pol}$, $\alpha_d=0.5$ (tesis) | — | 1345.5 |
| NIST | 557.1 | 1410.0 |

El $2107.7$ que hoy aparece en la columna HF+CI de esa tabla salió del $\zeta$ literal $0.00624$ Ha; el
correcto es $0.00375$ Ha, un factor $1.66$. **La tabla pone dos convenciones de $\zeta$ incompatibles
en el mismo renglón**, y toda la narrativa se construyó sobre ese artefacto.

Con $\zeta$ consistente el resultado **sobrevive y mejora**: HF+CI *subestima* (−11%), y $V_{pol}$
*sube* el desdoblamiento hasta −4.6%. El propio Cuadro `tab:estanio_ci` de la tesis ya lo demuestra:
$\zeta_{5p}$ crece con $\alpha_d$ (0.01001 a $\alpha_d=3.0$, 0.01073 a $\alpha_d=4.0$) frente al
$0.00859$ de $\alpha_d=0$. El mecanismo escrito en `:265` es el equivocado: $V_{pol} = -\alpha_d/2r^4$
es **atractivo**, contrae el orbital, y además $dV_{pol}/dr = +2\alpha_d/r^5$ entra directamente al
integrando de $\zeta$. Ambos efectos **suben** $\langle 1/r^3\rangle$.

**Qué hacer:** rehacer el barrido de $\alpha_d$ para Ge y Sn con el CI completo y $\zeta$ calculado
(`calibrate_vpol.jl`, `calibrate_tin_vpol.jl`), rellenar `file_vpol` en `generate_table.jl`, y
reescribir `:263-265`, el caption de `:248` y el párrafo de `:292`. La conclusión final es más fuerte
que la anterior, no más débil.

### B2. Canonizar los `.jld2` de C y Si
Reemplazar los originales por los `_f2fix` y borrar el sufijo, o dejarlos y actualizar todas las rutas.
Elegir uno; hoy conviven los dos y es una trampa. (~10 min)

### B3. Una regeneración completa desde archivos canónicos
Correr `generate_table.jl` una sola vez de principio a fin, con el manifiesto de procedencia, y pegar
el bloque de procedencia como comentario junto a cada tabla del `.tex`. Cuesta ~40 min, dominados por
el estaño. **Ningún número del capítulo 6 debe sobrevivir sin pasar por aquí.**

### B4. La columna del estaño de `tab:niveles_energia_pesados`
Está vacía porque no existe `tin_rohf_results_av_R30.0.jld2`. Decidir: o se genera, o los cuatro
elementos se reportan sobre orbitales $^3P$ y se corrige el comentario contradictorio de
`np2_toy_ci.jl`, que dice que el CI corre mejor sobre promedio de configuración mientras el código usa
los `_3P_`.

### B5. La fila de ionización del germanio
`resultados.tex:305` da 7.40 eV para "HF + CI (Valencia)" contra 7.33 de HF puro: 0.0026 Ha de
correlación, exactamente lo que recuperaba el CI roto. Recalcular con el CI bueno (que da 0.0071 Ha) y
reescribir el párrafo del 3.1% residual en `:311`.

### B6. El $F^0$ del germanio en `tab:integrales_slater`
El $0.30137782$ no se reproduce con ningún cálculo. Rastrearlo o recomputar la celda. El renglón del
germanio de esa tabla además es un copy-paste dígito a dígito del silicio.

## Definición de terminado

El código está congelado cuando B1–B6 están cerrados y `generate_table.jl` corre de principio a fin sin
avisos. A partir de ese punto **no se toca código**: cualquier discrepancia que aparezca durante la
redacción se documenta como limitación, no se persigue.

Estimación: B2 y B6 son de minutos. B1 es el trabajo real (dos barridos de $\alpha_d$ con SCF por
punto, que es justo para lo que se hizo C-DIIS). B3 es tiempo de máquina desatendido. B4 y B5 son
decisiones cortas más una corrida.

## Lo que va al capítulo de trabajo futuro, no al de resultados

- La razón $\zeta_{calc}/\zeta_{NIST}$ cae monótonamente con $Z$ (1.460, 0.952, 0.898, 0.870) y cruza
  1 entre C y Si. **No es un error sistemático único.** Apunta a dos omisiones distintas de la
  truncación de Breit-Pauli: en los ligeros faltan los términos de dos cuerpos (espín-otra-órbita y
  espín-espín), que reducen el desdoblamiento y pesan como $\sim 1/Z_{ef}$; en los pesados falta la
  contracción relativista del orbital de valencia, que lo aumenta y crece con $Z$. Ese cruce monótono
  es un resultado en sí mismo y hay que presentarlo como tal.
- CI multi-referencia con el $ns$ activo, para el $^1S$.
- Los términos de dos cuerpos de Breit.
