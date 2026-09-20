# Siguientes pasos: del código congelado al manuscrito

Escrito el 2026-09-20, con el código congelado en `1f22fc9`. Este archivo es el punto de entrada
para una sesión nueva: dice qué hay hecho, qué se puede tocar y en qué orden escribir. El detalle
de los números está en `REVISION-RESULTADOS-2026-09-13.md`, y el de la regeneración, en
`PLAN-REGENERACION.md`.

---

## 0. Dónde estamos

- **El cálculo está cerrado.** Manifiesto 25/25, 39 corridas de CI (carbono hasta m = 40), etapa 3
  reproducible y sin celdas faltantes. La regla de congelamiento está en `CLAUDE.md`.
- **Todo lo que la tesis reporta ya existe generado:** 16 fragmentos de tabla y `valores_texto.md`
  en `docs/tablas/`, y 12 figuras del pipeline en `docs/figures/`.
- **Lo que falta es el manuscrito.** Las ocho tablas del capítulo de resultados siguen escritas a
  mano, ninguna usa todavía `\input{tablas/...}`, y quedan 21 `\todo` en los `.tex`.

## 1. Qué hay listo para usar

| fragmento | dónde va |
|---|---|
| `energia_global.tex` | `tab:resultados_energia_global` |
| `momentos_inversos.tex` | `tab:momentos_inversos` |
| `integrales_slater.tex` | `tab:integrales_slater` |
| `koopmans.tex` | `tab:koopmans` |
| `niveles_ligeros.tex` | `tab:niveles_energia_ligeros` |
| `niveles_pesados.tex` | `tab:niveles_energia_pesados` |
| `barrido_vpol_Sn.tex` | `tab:estanio_ci` |
| `ionizacion_Ge.tex` | `tab:germanio_ionizacion` |
| `barrido_vpol_Ge.tex`, `ionizacion_Sn.tex`, `convergencia_ci.tex`, `convergencia_ci_costo.tex`, `convergencia_singletes.tex`, `cdiis.tex`, `zeta.tex`, `lande.tex` | secciones nuevas (punto 3) |

Los fragmentos son solo el `tabular`, con su procedencia en comentarios; la leyenda y la etiqueta
se quedan en el capítulo:

```latex
\begin{table}[htbp]
    \centering
    \input{tablas/energia_global}
    \caption{...}
    \label{tab:resultados_energia_global}
\end{table}
```

**Figuras generadas que el texto todavía no usa:** `convergencia_ci.pdf`,
`convergencia_singletes.pdf`, `cdiis_convergencia.pdf`, `zeta_razon.pdf` y `barrido_vpol.pdf`. Las
siete que ya citaba el capítulo se regeneraron con la misma base y los mismos nombres.

**`valores_texto.md`** tiene los números que la prosa cita fuera de las tablas: mínimos de V_rad,
ζ y F² en cm⁻¹, los ajustes de ζ al NIST, el CI por elemento, la sensibilidad de orbitales, la
convergencia de los singletes con sus extrapolaciones y el error de cada punto de los barridos.

**Diagnósticos disponibles** (solo lectura, segundos, en `benchmarks/np2_sequence/tesis/`):
`diagnostico_rayleigh.jl`, `diagnostico_casi_degeneracion.jl`, `diagnostico_malla_silicio.jl` y
`analisis_zeta.jl`.

## 2. Reglas mientras se redacta

- Vale el registro de redacción de `CLAUDE.md`: sinodales físicos, párrafos de al menos tres
  oraciones, nada de vocabulario inflado, toda afirmación fuerte con `\ref`, y comillas de Babel.
- El usuario compila el LaTeX; aquí solo verificación estática.
- **Ningún número se escribe a mano.** Si el texto necesita uno que no está, se agrega a
  `valores_texto.md` desde la etapa 3 y se cita de ahí.
- Una discrepancia nueva se documenta como limitación (punto 5), no se persigue con código.

## 3. Orden de trabajo

1. ~~**`resultados.tex`.**~~ **Hecho el 2026-09-20.** Las quince tablas entran por `\input`, no
   queda ningún `tabular` a mano, y los siete pasajes del punto 4 están reescritos. El capítulo
   pasó de 322 a 351 líneas.
2. ~~**Secciones nuevas en resultados.**~~ **Hecho el 2026-09-20.** Quedaron como
   `sec:cdiis_resultados` (dentro de la validación, que es donde pertenece un resultado del SCF),
   `sec:convergencia_ci`, `sec:zeta_diagnostico` y `sec:lande`. `convergencia_ci_costo` se fue al
   apéndice (`sec:costo_ci`), porque el registro de redacción manda el rendimiento ahí.
3. **`ciclo-scf.tex`** (6 `\todo`) y **`correlacion-estructura-fina.tex`** (1): el mecanismo de
   V_pol, la formulación de C-DIIS y el método de CI.
4. **`discusion.tex`:** reescritura completa. Son 73 líneas y cuatro de sus afirmaciones cambian.
5. **`conclusiones.tex`:** está vacío (7 líneas), y con él la sección de trabajo futuro, que ya
   tiene sus números en las secciones 6 y 7 de la revisión.
6. **Cierre:** introducción, prólogo, los 21 `\todo` restantes y una pasada de referencias.

## 4. Lo que los números contradicen

El detalle, con los pasajes exactos, está en la sección 3 de `REVISION-RESULTADOS-2026-09-13.md`.
En resumen:

| pasaje | dice hoy | dicen los números |
|---|---|---|
| `resultados.tex:248, 263-265` | HF sobrestima el desdoblamiento y V_pol lo baja | HF+CI lo subestima (−11 % Ge, −15 % Sn); V_pol sube ζ |
| `resultados.tex:292` | α_d = 3 acierta ³P₂ y 4 el ³P₁; compromiso en 3.5 | los cruces son 2.87 y 3.22; con α_d = 3.0 los dos quedan a menos de 2 % |
| `resultados.tex:294-311` | V_pol reduce a la mitad el error de ionización | lo empeora: 1.4 % → 4.8 % en germanio |
| `resultados.tex:271, 311` | el ¹S₀ alto es del espacio CI truncado | el espacio está convergido; es del modelo de pareja con core congelado |
| `discusion.tex:64-73` | g = 1.4975 frente a 1.496, y el CI cierra la brecha | el ASD da 1.49458: el germanio mezcla 1.31 % de ¹D₂ y el modelo 0.8-0.9 % |
| ~~`resultados.tex:172`~~ | la figura confirma ζ ∝ Z⁴ | resuelto: el ajuste ya sale de la etapa 3 (p = 1.82 contra el Z nuclear) y el pasaje explica que la ley Z_eff⁴ es sobre la carga efectiva |
| `ciclo-scf.tex:205` | NIST ASD versión 6.1 | es la 5.12; citar con la clave `NIST_ASD` de la bibliografía |

Además, las tablas ya corrigieron el virial del estaño, el ⟨r⁻³⟩ del germanio y el F⁰(4p,4p) del
germanio: el texto que los comenta tiene que seguir a las tablas, no al revés.

## 5. Limitaciones que el texto debe declarar

- **ζ de un cuerpo.** En los ligeros faltan los términos de dos cuerpos de Breit-Pauli; en los
  pesados el NIST pide más que la cota de núcleo desnudo, así que falta contracción del orbital,
  relativista y/o por polarización del core (`PLAN-REGENERACION.md` §6).
- **Singletes.** El CI de pareja con core congelado no contiene la casi-degeneración ns² ↔ np²; el
  tamaño del efecto y su no aditividad están medidos (revisión §7).
- **α_d efectivo.** Calibrado contra el NIST absorbe algo más que polarización: empeora la
  ionización (revisión §2.5).
- **Ortogonalización por Gram-Schmidt** en vez de multiplicadores de Lagrange: explica el 0.2-0.3 %
  de ⟨r⁻³⟩ y la desviación del virial frente a Froese Fischer (revisión §8).
- **Energía con autovalores** en vez del cociente de Rayleigh (`PLAN-REGENERACION.md` §2.2).
- **Truncamiento del CI en lmax = 3.**

## 6. Decisiones abiertas

- ~~**Modo de trabajo.**~~ **Decidido el 2026-09-20:** Claude redacta los `.tex` y el usuario
  revisa, capítulo por capítulo en el orden del punto 3.
- **ζ(2p) = 31.946 cm⁻¹ de Froese Fischer:** no se sabe si es el de Blume-Watson. Cambia cómo se
  cita en la discusión de ζ.
- **Polarizabilidades de los cores de Ge y Sn:** sin fuente verificada (revisión §6.3). Sin ellas,
  la comparación con α_d queda cualitativa.
- **Casi-degeneración:** queda como trabajo futuro (opción A) salvo que el usuario decida invertir
  semanas en un CI de cuatro electrones.

## 7. Cómo arrancar una sesión nueva

```bash
julia benchmarks/verificar_resultados.jl
julia --project=. benchmarks/np2_sequence/tesis/etapa2_ci.jl --lista
```

El primero tiene que decir 25/25 y el segundo, todo `listo`. Si alguno se queja, algo se movió y
hay que entenderlo antes de escribir una sola línea. Después: este archivo, las secciones 3 y 5 de
`REVISION-RESULTADOS-2026-09-13.md` y `docs/tablas/valores_texto.md`.
