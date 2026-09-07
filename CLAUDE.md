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

## Cómo correr

```bash
julia --project=. benchmarks/np2_sequence/np2_toy_ci.jl      # CI de pares (legado, verificado)
julia --project=. benchmarks/np2_sequence/np2_ci_full.jl     # CI general (sin ejecutar aún)
julia --project=. benchmarks/np2_sequence/tin_rohf.jl        # SCF del estaño (~200 iteraciones)
```

Los `.jld2` de geometría se cachean por `(Z, R_max, N, K, γ, alpha_d, r_c)`; borrar solo si cambia
la malla. Los `*_rohf_results_*.jld2` guardan `orbitals`, `E_total`, `R_grid`, `V_eff`, `P_np`.

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
  closed_shell/         átomos de capa cerrada; 11_radon.jl tiene el prototipo de C-DIIS que funciona
  np2_sequence/         la secuencia C, Si, Ge, Sn de la tesis
docs/
  main.tex, chapters/, appendices/, claude/   manuscrito y notas de trabajo
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
