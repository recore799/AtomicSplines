# `benchmarks/` — qué se versiona aquí y por qué

Este directorio mezcla scripts de cálculo con los datos que producen. Los datos
no son todos iguales y no se tratan igual.

## Las tres clases de archivo

**1. Cachés de geometría (`geometry_*.jld2`, de 0.9 a 16 MB).** Puramente
derivadas: se regeneran desde `assemble_geometry` con la clave
`(Z, R_max, N, K, γ, alpha_d, r_c)`. Viven en `.geometry_cache/` en la raíz del
repo y están ignoradas por git. Borrarlas solo cuesta tiempo de máquina. Si
cambia la malla, se borran y ya.

**2. Resultados de SCF (`*_rohf_results_*.jld2`, de 0.3 a 2.3 MB).** Guardan
`orbitals`, `E_total`, `R_grid`, `V_eff` y `P_np`. **Están versionados en git**
y registrados en [`RESULTS.toml`](RESULTS.toml). Son la entrada de las tablas
del capítulo de resultados: regenerar uno cuesta desde segundos (C, Si) hasta
diez o quince minutos (Sn), pero el problema nunca fue el costo sino no poder
fechar un cambio.

**3. Scratch.** No existe. Si vuelve a aparecer, o tiene un nombre que dice para
qué sirve, o se borra.

## Por qué los resultados sí van a git

El 2026-09-06 se encontró que `carbon_rohf.jl` y `silicon_rohf.jl` tenían
invertido el signo del acoplamiento `coeff_k2` del término ³P, y que el código
había regresionado *después* de generar las tablas del capítulo 6. Los `.jld2`
de resultado no estaban trackeados, así que no hubo `git log -p` que señalara el
commit culpable: la evidencia hubo que reconstruirla comparando a mano contra
Froese Fischer y contra los números ya impresos en `resultados.tex`. El detalle
está en `docs/claude/HALLAZGOS-2026-09-06.md`.

De las tres opciones sobre la mesa —git directo, `git-lfs`, o solo un
manifiesto— se tomaron dos, y esta es la razón:

- **Git directo, sin LFS.** El conjunto completo pesa ~14 MB. Git lo maneja sin
  problema a este tamaño, y meter `git-lfs` obligaría a instalarlo a cualquiera
  que clone el repo solo para poder leer la entrada de las tablas de la tesis,
  que es exactamente lo contrario de lo que se busca. Si el conjunto creciera un
  orden de magnitud, la decisión se revisa.

- **Manifiesto encima.** Un `.jld2` versionado aparece en `git status` cuando
  cambia, pero `git diff` sobre un binario solo dice «Binary files differ».
  `RESULTS.toml` fija por archivo el sha256, la energía total, los parámetros
  físicos `(Z, R_max, N, K, γ, alpha_d, r_c)`, el estado (`av` / `3P`), el
  script que lo generó y el commit de ese script. Así el diff **se lee**:

  ```
  -E_total       = -37.594380215677
  +E_total       = -37.6886189635015
  ```

El criterio de éxito es concreto: si alguien vuelve a invertir ese signo, no
hace falta que nadie se acuerde de comparar contra Froese Fischer. El resultado
cambia, el sha256 deja de cuadrar y la verificación falla sola.

## Verificar

```bash
julia benchmarks/verificar_resultados.jl
```

Compara el sha256 de cada resultado contra el manifiesto y avisa de archivos que
falten o que nadie haya registrado. No carga el proyecto —solo usa `SHA` y
`TOML` de la biblioteca estándar—, así que sirve igual en un hook de `pre-commit`
o en CI.

Si regeneraste un SCF a propósito, reemite el manifiesto y commitea el diff:

```bash
julia --project=. benchmarks/verificar_resultados.jl --emitir > benchmarks/RESULTS.toml
```

Los metadatos (parámetros, notas, script generador) viven en el `CATALOGO` de
`verificar_resultados.jl`, no en el `.toml`; reemitir no los pierde. Un resultado
nuevo se registra agregándolo a ese catálogo.

## Dónde corren los scripts

Los `*_rohf.jl` escriben su `.jld2` con ruta relativa **al directorio del
script**, no al directorio de trabajo, y los consumidores leen del mismo lugar.
Se puede lanzar desde donde sea:

```bash
julia --project=. benchmarks/np2_sequence/tin_rohf.jl
```

Antes no era así, y por eso había resultados sueltos en la raíz del repo, entre
ellos los cinco del barrido de $\alpha_d$ del estaño que sostienen el Cuadro
`tab:estanio_ci`.

## Nomenclatura

`<elemento>_rohf_results_<estado>_R<R_max>[_ad<alpha_d>].jld2`, con `estado` en
`{av, 3P}`. El sufijo `_ad` solo aparece cuando hay potencial de polarización del
core. **Un archivo por elemento y estado**: no se acumulan sufijos de parche
(`_f2fix` y compañía). Si un resultado se corrige, se reemplaza el archivo y el
historial de git guarda el anterior.

## Estado abierto

- **Falta `tin_rohf_results_av_R30.0.jld2`** (estaño, promedio de configuración,
  sin `V_pol`). Su ausencia es la razón de que la columna del estaño de
  `tab:niveles_energia_pesados` esté vacía, y de que `tin_validate.jl` caiga
  siempre a la rama del ³P. No es intencional: hay que generarlo.

  ```bash
  julia --project=. -e 'include("benchmarks/np2_sequence/tin_rohf.jl"); solve_tin_rohf(30.0; estado="av", use_diis=true)'
  ```

  Son ~10-15 minutos. Después, reemitir el manifiesto.

- **`carbon_rohf_results_3P_R30.0.jld2` y su equivalente de silicio** se
  escribieron a mano en la sesión de la corrección y les faltan las claves
  `active_s` y `active_p` que sí guarda el `jldsave` de los scripts. Ningún
  consumidor de `np2_sequence/` las lee, pero el esquema quedó incompleto.
  Regenerarlos con el script actual cuesta segundos y lo cierra.

- **`np2_toy_ci.jl` corre sobre los archivos `_3P_`** mientras su propio
  comentario dice que el CI se hace mejor sobre los de promedio de
  configuración. El comentario y el código se contradicen; hay que decidir cuál
  gana.
