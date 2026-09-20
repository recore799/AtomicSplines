# AtomicSplines.jl — contexto permanente

Librería en Julia para estructura atómica con base de B-splines. Es el motor de cálculo de la
tesis de licenciatura del usuario; el manuscrito vive en `docs/`.

## Reglas duras

1. **El usuario compila el LaTeX.** Nunca ejecutar `latexmk`, `pdflatex` ni `bibtex`. Verificación
   estática únicamente (leer, hacer grep, comprobar `\label`/`\ref`).
2. **Nunca cambiar un número de una tabla de la tesis sin regenerarlo y dejar registro de su
   procedencia.** Si un número no se puede regenerar con el código actual, decirlo en vez de copiarlo.
3. **Nunca hardcodear una cantidad física que el código puede calcular.** Este repo ya sufrió por eso:
   había un $\zeta$ literal en un sitio y calculado en otro, y tres valores distintos de $F^2$ para el
   mismo átomo repartidos entre un script, una tabla y una constante.
4. Los cálculos pesados (SCF de Sn, CI grandes) tardan minutos. Avisar antes de lanzarlos y guardar
   los logs en archivo.
5. **Los `*_rohf_results_*.jld2` están versionados y registrados en `benchmarks/RESULTS.toml`.**
   Si se regenera uno, hay que reemitir el manifiesto y commitear el diff; si no, se rompe la
   única señal que detecta una regresión numérica sola. Ver `benchmarks/README.md`.

## Congelamiento del código (desde 2026-09-13)

El código que produce los números de la tesis está congelado: manifiesto 25/25, 39 corridas de CI
registradas y etapa 3 reproducible. La fase en curso es la redacción del manuscrito, y el orden de
trabajo está en `docs/claude/SIGUIENTES-PASOS.md`.

**Se puede sin preguntar:**

- editar los valores de literatura de `benchmarks/np2_sequence/tesis/config.jl` y regenerar la
  etapa 3 (segundos);
- agregar a `valores_texto.md` una salida nueva que el texto necesite citar;
- escribir diagnósticos de solo lectura, con nombre propio, que no toquen los resultados.

**No se toca sin que el usuario lo pida explícitamente:** la física del SCF, del CI o de
Breit-Pauli; los `.jld2` y el `resultados_ci.toml`; y las decisiones de `config.jl` (`ESTADO`,
`ENERGIA_HF`, `CRITERIO_ALPHA`, `TAMANOS_CONVERGENCIA`, `BARRIDO_ALPHA`).

Una discrepancia que aparezca mientras se redacta se documenta como limitación; no se persigue con
código. Las que ya están diagnosticadas viven en `docs/claude/REVISION-RESULTADOS-2026-09-13.md`.

## Cómo correr

```bash
julia benchmarks/verificar_resultados.jl                              # los .jld2 siguen siendo los registrados
julia --project=. benchmarks/np2_sequence/tesis/etapa1_scf.jl         # qué SCF faltan (--correr los corre)
julia --project=. benchmarks/np2_sequence/tesis/etapa2_ci.jl --lista  # qué CI faltan (sin --lista, ~2.5 h)
julia --project=. benchmarks/np2_sequence/tesis/etapa3_tablas.jl      # docs/tablas/, segundos
julia --project=. benchmarks/np2_sequence/tesis/etapa4_figuras.jl     # docs/figures/, ~1 min
julia --project=. benchmarks/np2_sequence/np2_ci_full.jl              # CI general, prueba sobre carbono
```

Las decisiones de la regeneración (orbitales, `alpha_d`, espacio activo, fórmula de la energía) y
los datos de literatura (NIST, Froese Fischer) viven solo en `benchmarks/np2_sequence/tesis/config.jl`.
El orden completo, los tiempos y qué revisar están en `docs/claude/PLAN-REGENERACION.md`.

Las rutas de datos están ancladas al directorio del script, no al de trabajo: da igual desde dónde
se lance `julia`.

Los `.jld2` de geometría se cachean por `(Z, R_max, N, K, γ, alpha_d, r_c)` en `.geometry_cache/`
en la raíz, ignorada por git; borrar solo si cambia la malla. Los `*_rohf_results_*.jld2` guardan
`orbitals`, `E_total`, `R_grid`, `V_eff`, `P_np` (y `alpha_d`, `r_c` cuando llevan V_pol), y **sí**
se versionan: son la entrada de las tablas del capítulo 6.

## Mapa

```
src/
  structs.jl            Orbital{n,l,occ,energy,coeffs}, SolverWorkspace
  basis.jl              B-splines, malla exponencial
  integrals.jl          assemble_geometry (incluye V_pol en ws.V), Fock J/K, cached_init_scf_workspace
  poisson.jl            solve_generalized_poisson!
  atomic_calculations.jl compute_effective_central_potential
  ci.jl                 compute_Rk, get_cached_Rk!, extract_virtuals, CFPs
benchmarks/
  README.md             qué se versiona aquí y por qué; nomenclatura de los datos
  RESULTS.toml          manifiesto: sha256, E_total, parámetros y procedencia de cada resultado
  verificar_resultados.jl  compara los .jld2 contra el manifiesto (--emitir lo regenera)
  closed_shell/         átomos de capa cerrada; 11_radon.jl tiene el prototipo de C-DIIS que funciona
  np2_sequence/         la secuencia C, Si, Ge, Sn de la tesis
                        <el>_rohf.jl producen los *_rohf_results_*.jld2; Ge y Sn aceptan alpha_d
                        np2_ci_full.jl es el motor de CI (T1-T4 se autoverifican)
                        resultados_ci.toml: los CI que alimentan las tablas (lo escribe la etapa 2)
                        tesis/  config.jl y las cuatro etapas de la regeneración
                        legado, ya sin números en la tesis: el CI de pares de np2_toy_ci.jl (sus
                        utilidades, como compute_zeta, las usa np2_ci_full.jl), <el>_rohf_vpol.jl,
                        generate_table.jl, calibrate_*.jl, *_fine_structure.jl, scaling_law.jl,
                        plot_energy_levels.jl y examples/scratch/plot_*.jl
.geometry_cache/        cachés de geometría, derivadas e ignoradas por git
docs/
  main.tex, chapters/, appendices/, claude/   manuscrito y notas de trabajo
  tablas/               fragmentos tabular de la etapa 3; no se editan a mano
```

## Registro de redacción (para cualquier cosa que toque `docs/*.tex`)

- Los sinodales son físicos. Maquinaria teórica en el mainmatter; rendimiento, BLAS, caché y memoria
  en el apéndice. **El capítulo muestra los resultados; el apéndice, los pasos.**
- Párrafos completos de al menos tres oraciones. Nunca oraciones sueltas.
- Prohibido el vocabulario inflado ("triunfo fundamental", "robustez excepcional", "masivo",
  "vertiginoso"): lee como redacción de IA.
- Toda afirmación fuerte cierra con una referencia interna (`\ref`). Nada de contenido duplicado:
  cada ecuación se define en un solo lugar.
- Babel va con `es-noquoting`: nunca comillas rectas, usar ``...''.
