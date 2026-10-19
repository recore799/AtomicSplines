# Auditoría final del manuscrito, en seco

Haz una auditoría completa del manuscrito de la tesis antes de enviarlo al asesor. Es **en seco**:
no modifiques ningún `.tex`, `.jl`, `.toml`, `.bib` ni `.jld2`, no regeneres tablas ni figuras y no
hagas commits. El único archivo que escribes es el reporte. Trabaja tú solo, sin subagentes, porque
la mitad de los hallazgos son de coherencia entre capítulos y requieren una sola lectura del todo.

## Qué se audita

Lo que incluye `docs/main.tex`, en su orden: resumen, prólogo, capítulos 1 a 7 y apéndices A a E,
incluidos los fragmentos de `docs/tablas/` que se insertan con `\input` y los captions. Quedan
fuera `docs/chapters/_borradores/` y las tablas `*_vpol.tex`, que la tesis ya no usa.

## Cómo trabajar

1. **Primero lee todo de corrido**, como lo leería un sinodal: un físico que conoce Hartree-Fock al
   nivel del Szabo y trabaja en sistemas confinados con DFT. Toma notas, pero no reportes nada hasta
   terminar la lectura.
2. **Después haz las revisiones dirigidas** de la lista siguiente, con `grep` cuando se pueda contar.
3. **Compila en un directorio temporal**, sin ensuciar `docs/`. Desde `docs/`, ejecuta
   `latexmk -pdf -interaction=nonstopmode -file-line-error -outdir=<scratch> main.tex` y guarda el
   log en archivo. Mira como imagen (`pdftoppm`) las páginas con tablas anchas o figuras grandes.

## Qué buscar, por gravedad

**[A] Errores.** Un número de la prosa que no coincide con la tabla generada que lo respalda
(`docs/tablas/*.tex` o `docs/tablas/valores_texto.md`). Una afirmación que contradice otra parte del
texto. Un signo, un factor o un índice equivocado en una ecuación. Una afirmación sobre lo que hace
el código que no se sostiene al leer `src/` o `benchmarks/np2_sequence/`; verifícalas leyendo el
código, no las notas de `docs/claude/`. Una referencia que apunta a la sección equivocada.

**[B] Lo que un sinodal puede objetar.** Afirmaciones fuertes sin `\ref` interna ni cita, o con
más alcance del que dan los números. Términos usados antes de definirse. Saltos de razonamiento.
Atribuciones dudosas a una fuente. Las citas con localizador (sección, capítulo, página) no las
puedes comprobar sin el libro, así que no las des por buenas ni por malas: lístalas al final para
que el autor las coteje.

**[C] Redacción.** Hay que quitar los patrones que delatan texto generado por un LLM sin imponer
otro estilo, porque el autor quiere que suene a él:
- «no es X sino Y», «no solo… sino»;
- tríadas de adjetivos, de verbos o de ejemplos;
- cadenas de punto y coma y de dos puntos;
- frases finales que resumen el párrafo que acaban de cerrar;
- conectores repetidos («de modo que», «es decir», «en cambio», «así», «por tanto», «precisamente»);
- arranques anafóricos en serie («Ese…», «Esa…», «Esto…»);
- párrafos de largo uniforme y vocabulario inflado.

Además, revisa los calcos del inglés (como lo era «carácter ¹D₂») y la ortografía: tildes, «solo» y
los demostrativos sin tilde, y concordancia. Respeta el registro de `CLAUDE.md`:
- párrafos de al menos tres oraciones;
- nada de oraciones sueltas;
- cada ecuación definida en un solo lugar;
- maquinaria teórica en el cuerpo y rendimiento en el apéndice;
- comillas ``...'' (con `es-noquoting`).

Para cada caso propone la reescritura mínima que resuelve el problema y no cambia números,
referencias, ecuaciones ni afirmaciones físicas. No reportes cada aparición de un tic: da el conteo
global y los casos en que de verdad estorba.

**[D] Tipografía y formato.** Cajas desbordadas, referencias o citas indefinidas, flotantes lejos
de su primera mención o que dejan páginas medio vacías, espacios entre número y unidad, `\%` y
`\text` en modo matemático, coherencia de notación entre capítulos (ζ_np, F^k, ⟨r⁻³⟩, E_h, cm⁻¹).
Las observaciones que afectan a todo el documento, como el uso de mayúsculas en los títulos, van
una sola vez en la sección de observaciones globales.

## Decisiones ya tomadas: no las reabras

- El potencial de polarización (V_pol) está fuera del modelo. Sus menciones en el final de la
  Sección 6.2, en la 6.6.3 y en la 7.2 son intencionales.
- Terminología: «capas internas» (nunca «core»); «triplete»; «componente ¹D₂»; «núcleo» solo para
  el núcleo atómico; «núcleo desnudo», definido en la Sección 3.2.
- ζ_np se calcula; no se ajusta. Que el código multiplique los elementos de espín-órbita por los
  coeficientes CI ya está documentado en el Capítulo 5.
- Las discrepancias con el experimento ya diagnosticadas
  (`docs/claude/REVISION-RESULTADOS-2026-09-13.md` y la Sección 6.6) son limitaciones documentadas.
  No se persiguen con código.
- Formato propio a 11 pt, «Tabla» en lugar de «Cuadro», orden actual de las secciones del Capítulo 6.
- «de la de la» en el paso 4 del ciclo SCF es intencional.

Si algo de esto te parece un error de verdad, dilo una vez en las observaciones globales, con el
argumento. No lo repitas hallazgo por hallazgo.

## Entregable

Escribe `docs/revision/auditoria-final.org` en el formato de `docs/revision/observaciones.org`, que
es el archivo donde el autor marca sus pendientes. Debe tener:

1. **Veredicto** en un párrafo y los diez hallazgos que corregirías antes de enviar.
2. **Observaciones globales**, incluido el conteo de cada tic de redacción.
3. **Hallazgos por capítulo y apéndice**, cada uno como `** TODO [A|B|C|D] §x.y — resumen`, con:
   - `archivo:línea`;
   - la cita textual mínima;
   - el problema en una oración;
   - la propuesta.
4. **Citas con localizador** que el autor debe cotejar con el libro.
5. **Resultado de la compilación** (páginas, advertencias) y las páginas que miraste como imagen.

Todo en español. Al terminar, resume en el chat cuántos hallazgos hay de cada gravedad y cuáles son
los tres más importantes. No corrijas nada hasta que el autor diga cuáles aplicar.
