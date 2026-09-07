# Prompt: organizar y trackear los resultados de benchmarks/

Lee `CLAUDE.md` completo antes de tocar nada — tiene reglas duras de este repo
(nunca cambiar un número de tabla de la tesis sin regenerarlo con procedencia,
nunca hardcodear una cantidad física que el código calcula).

## Por qué esto importa ahora

`benchmarks/np2_sequence/` mezcla tres cosas con reglas de vida distintas, y hoy
todas viven igual de desprotegidas porque `.gitignore` tiene `*.jld2` a nivel
global:

1. **Cachés de geometría puramente derivadas** (`geometry_Z*.jld2`, de
   6 a 16 MB cada una): se regeneran desde `assemble_geometry` a partir de la
   clave `(Z, R_max, N, K, γ, alpha_d, r_c)`. Estar sin trackear es correcto.
2. **Resultados de SCF que alimentan la tesis** (`*_rohf_results_*.jld2`,
   300 KB–2.3 MB): `orbitals`, `E_total`, `V_eff`. Estos son los que de verdad
   importan y hoy están completamente indefensos.
3. **Basura de sesiones pasadas**: `.bak`, `*_output.txt`, scripts de prueba
   sueltos (`test_ci.jl`, `test_ci2.jl`, `fix_script.py`, etc.), duplicados de
   geometría con y sin sufijo `_ad0.000_rc1.000.jld2` para los mismos parámetros.

**El incidente que expuso el problema**: el 2026-09-06 se encontró que
`carbon_rohf.jl` y `silicon_rohf.jl` tenían un signo de acoplamiento invertido
(`coeff_k2` del término ³P). El código regresionó *después* de que se generaron
las tablas del capítulo 6 — se pudo probar comparando contra Froese Fischer y
contra los propios números ya impresos en `resultados.tex`, no contra ningún
historial del repo, porque **los `.jld2` de resultado nunca estuvieron
trackeados**. Si lo hubieran estado, `git blame`/`git log -p` en el `.jld2`
habría mostrado el commit exacto que rompió el signo, con fecha y diff. En vez
de eso hubo que reconstruir la evidencia por triangulación numérica. El detalle
completo está en `docs/claude/HALLAZGOS-2026-09-06.md` — léelo, tiene el
inventario de qué se corrigió y qué archivos quedaron contaminados.

Ahora mismo hay dos archivos de resultado con dos versiones alternativas
conviviendo, `carbon_rohf_results_3P_R30.0.jld2` (contaminado por el bug) junto
a `carbon_rohf_results_3P_R30.0_f2fix.jld2` (correcto), y lo mismo para
silicio. Esa situación tiene que resolverse como parte de este trabajo, no
quedar así.

## Lo que hay que decidir y hacer

1. **Clasificar cada archivo no-`.jl` de `benchmarks/np2_sequence/`** (y de
   `benchmarks/closed_shell/` si aplica el mismo patrón) en una de las tres
   categorías de arriba. Usa `ls -la` para ver tamaños reales antes de asumir.

2. **Elegir cómo trackear los resultados de SCF que alimentan la tesis**
   (categoría 2) y justificar la elección en un comentario o README corto
   dentro de `benchmarks/`:
   - Trackearlos directamente en git (son pocos MB, git los maneja bien sin
     LFS a este tamaño), o
   - Usar `git-lfs` si prefieres no meter binarios al historial normal, o
   - Como mínimo, un manifiesto versionado (JSON/TOML) por archivo con hash
     sha256, el commit del código que lo generó, los parámetros físicos
     (Z, R_max, N, K, γ, alpha_d, r_c, estado 3P/av) y la energía total, de
     modo que un `git diff` sobre el manifiesto detecte una regresión aunque
     el binario en sí no esté versionado.

   Cualquiera que elijas, el criterio de éxito es: **si alguien vuelve a
   invertir ese signo por accidente, `git status`/`git diff`/`git log` tiene
   que mostrarlo**, no requerir que alguien recuerde comparar contra Froese
   Fischer a mano.

3. **Resolver los archivos `_f2fix`**: decide si reemplazan a los originales
   contaminados (mi recomendación en `HALLAZGOS-2026-09-06.md` es que sí, una
   vez que confirmes el razonamiento del signo) y deja un solo archivo por
   elemento/estado, sin sufijos de parche acumulándose.

4. **Resolver la ambigüedad de nombres**: existen `carbon_rohf_results_R30.0.jld2`
   (sin sufijo de estado), `..._av_R30.0.jld2` y `..._3P_R30.0.jld2` para el
   mismo elemento. Averigua si el primero es legado muerto de antes de que el
   ciclo se dividiera en promedio de configuración / estado específico, y si lo
   es, elimínalo. Nota también que **no existe** `tin_rohf_results_av_R30.0.jld2`
   — decide si hace falta generarlo o si la ausencia es intencional.

5. **Deduplicar las cachés de geometría**: hay pares como
   `geometry_Z32.0_R30.0_N300_K8_g3.00.jld2` y
   `geometry_Z32.0_R30.0_N300_K8_g3.00_ad0.000_rc1.000.jld2` que probablemente
   son la misma malla con dos convenciones de nombre distintas (antes y
   después de que se agregara `alpha_d`/`r_c` a la clave del caché). Verifica
   si son idénticas (mismo hash o mismas dimensiones/valores) antes de borrar
   ninguna — no asumas.

6. **Limpiar basura de sesiones pasadas**: `*.bak` sueltos en el working tree
   (ya cubiertos por `.gitignore` pero ensuciando `ls`), archivos `*_output.txt`
   con volcados de stdout, y los scripts `test_*.jl`/`fix_script.py` de una
   sola sesión. Para cada uno: o tiene un propósito claro y permanente (dale un
   nombre que lo diga y consérvalo), o es scratch de una sesión pasada
   (bórralo). No lo dejes en ambigüedad de nuevo.

7. **Actualizar `.gitignore`** para reflejar la decisión del punto 2: si los
   resultados de SCF pasan a trackearse, la regla global `*.jld2` tiene que
   volverse más específica (por ejemplo, ignorar solo `geometry_*.jld2` o
   moverlas a un subdirectorio tipo `.cache/` que sí esté ignorado por completo).

## Cómo trabajar

- Este es un cambio de infraestructura, no de física: **no toques ningún
  número de las tablas de la tesis ni el contenido de ningún `.jld2` bueno**,
  solo su organización y su estado de tracking.
- Trabaja en commits pequeños y explicables (uno por decisión: "trackear
  resultados SCF", "eliminar cachés de geometría duplicadas", "limpiar scripts
  de prueba muertos"), no un solo commit masivo.
- Antes de borrar cualquier `.jld2`, confirma que no es la única copia de algo
  irremplazable — corriendo el SCF correspondiente de nuevo cuesta desde
  segundos (C, Si) hasta ~10-15 minutos (Sn), así que no es catastrófico
  regenerarlo, pero avisa antes de un `rm` sobre un resultado que no puedas
  reproducir en el momento.
- Al final, entrega un resumen corto de qué se trackeó, qué se eliminó, qué
  se dejó igual y por qué, y actualiza el mapa del repo en `CLAUDE.md` si la
  estructura de `benchmarks/np2_sequence/` cambió de forma relevante para
  quien trabaje aquí después.
