# Revisión de resultados — 2026-09-13

Revisión de las cuatro etapas que corrió el usuario el 13-09 (commits `d267ed1`, `b74f25a` y
`af9e882`), siguiendo `PLAN-REGENERACION.md` §4. Los números salen de
`benchmarks/np2_sequence/resultados_ci.toml` y `docs/tablas/`. Los derivados (errores, cruces,
extrapolaciones) dicen cómo se obtuvieron.

---

## 0. En corto

- **Las corridas están sanas y el código se puede congelar.**
  - Los 9 SCF convergieron con C-DIIS y las 18 autoverificaciones del CI pasan.
  - A m = 20 el HF+CI reproduce exactamente el log del 07-09.
  - La etapa 3 es reproducible y solo deja vacías las celdas de literatura.
- **Queda un punto abierto del lado del cálculo:** los singletes del carbono no están
  convergidos a m = 20 (sección 2.2).
- **Cuatro afirmaciones de la tesis se invierten con los números nuevos** (sección 3): el sentido
  de V_pol, la ionización con V_pol, la brecha del factor g del germanio y la tensión del barrido
  del estaño.

---

## 1. Controles

| control | resultado |
|---|---|
| SCF de la etapa 1 | 9/9 con C-DIIS, sin caer al respaldo; Ge en 65–69 iteraciones (~1.5 min), Sn en 35–39 (~3.5 min) |
| manifiesto | 25/25; el diff solo agrega las 9 entradas |
| autoverificación del CI | 18/18 pasan T1–T4 |
| raíces | 34/34 en [1, 1, 1] |
| procedencia | las 34 corridas, con el código de `2479092` sin cambios |
| reproducibilidad | a m = 20 y sin V_pol, E_corr y ³P₂ idénticos al log del 07-09 en los cuatro elementos |
| etapa 3 | sin celdas faltantes salvo las de literatura; volver a correrla da archivos idénticos |

---

## 2. Qué dicen los números

### 2.1 Estructura fina sin V_pol (m = 20, error frente al NIST)

| | ³P₁ | ³P₂ | ¹D₂ | ¹S₀ |
|---|---|---|---|---|
| C | +28.9 % | +45.3 % | +9.4 % | +13.5 % |
| Si | −6.2 % | −5.7 % | +15.7 % | −4.7 % |
| Ge | −11.9 % | −10.8 % | +10.2 % | −3.6 % |
| Sn | −19.2 % | −14.8 % | +1.1 % | −6.5 % |

Los ³P_J están convergidos en el espacio activo: el ³P₂ del carbono se mueve 0.13 cm⁻¹ de m = 4
a m = 20. Su error sigue la razón de ζ (`PLAN-REGENERACION.md` §6).

### 2.2 Convergencia de los singletes

Extrapolación geométrica con los incrementos de m = 12 → 16 → 20:

| | ¹D₂ a m = 20 | cola | ¹S₀ a m = 20 | cola |
|---|---|---|---|---|
| C | 11148.4 | −671 | 24574.8 | −870 |
| Si | 7285.2 | −31 | 14674.5 | −103 |
| Ge | 7851.6 | −20 | 15780.2 | −32 |
| Sn | 8710.1 | no fiable (razón 0.93) | 16045.3 | −45 |

- **Si, Ge y Sn están convergidos:** les falta del orden de 100 cm⁻¹ o menos.
- **El carbono no.** A m = 20 el ¹D₂ todavía baja 200 cm⁻¹ por paso y el ¹S₀, 530.
  - Extrapolados, ¹D₂ ≈ 10 480 cm⁻¹ (+2.8 %) y ¹S₀ ≈ 23 700 (+9.5 %).
  - Es el patrón que predice la casi-degeneración ns² ↔ np²: el ¹D₂ se acerca al NIST y el ¹S₀ se
    queda alto.
  - Pero es una extrapolación, y la tesis quiere afirmarlo (ver sección 4.1).

### 2.3 Sensibilidad ³P frente a promedio de configuración (m = 20)

Con orbitales de promedio, los ³P_J bajan 1.5–1.8 % (la razón de ζ) y los singletes 0.1–0.6 %.
La elección de orbitales no cambia ninguna conclusión.

### 2.4 V_pol

| | cruce ³P₁ | cruce ³P₂ | α_d elegido | ³P₁ | ³P₂ | ¹D₂ | ¹S₀ |
|---|---|---|---|---|---|---|---|
| Ge | 0.67 | 0.66 | 0.75 | +1.6 % | +1.5 % | +15.3 % | +0.6 % |
| Sn | 3.22 | 2.87 | 3.0 | −1.9 % | +0.8 % | +12.3 % | +2.1 % |

«Cruce» es el α_d en que el intervalo coincide con el NIST, interpolando linealmente el barrido;
los errores son con el α_d elegido.

- **V_pol sube ζ:** +13.5 % en Ge (α_d = 0.75) y +18.1 % en Sn (α_d = 3.0). F² sube menos, +2.6 %
  y +5.6 %: ζ es de 3 a 5 veces más sensible.
- **Germanio:** un solo α_d (≈ 0.67) ajusta los dos intervalos a la vez.
- **Estaño:** la tensión entre intervalos existe pero es chica: los cruces son 2.87 y 3.22, no 3 y
  4 como dice la tesis.
- **Los singletes sí se mueven.** El ¹D₂ sube 364 cm⁻¹ (Ge) y 965 (Sn), y el ¹S₀, 683 y 1485.
  - Parte del corrimiento es la referencia ³P₀, que baja al crecer ζ. La separación LS ¹D − ³P del
    CI crece 3.2 % y 6.6 %.
  - El ¹D₂ empeora y el ¹S₀ mejora.
- **Hipótesis para la discusión.** El ¹D₂ alto (ya sin V_pol en Si y Ge) apunta a un F² efectivo
  demasiado grande. Es lo que corregiría el término de dos cuerpos de la polarización del core,
  que el modelo omite (Apéndice F, ecuación `eq:cpp_dos_cuerpos`). Cowan (1981) reduce los F^k de
  Hartree-Fock en sus ajustes por efectos de correlación; verificar el factor y su justificación
  antes de citarlo.

### 2.5 Ionización (eV, orbitales ³P)

| | −ε | error | CI (ΔE) | error | + V_pol | error |
|---|---|---|---|---|---|---|
| C | 11.79 | +4.7 % | 12.03 | +6.8 % | | |
| Si | 8.08 | −0.8 % | 8.32 | +2.1 % | | |
| Ge | 7.82 | −1.0 % | 8.01 | +1.4 % | 8.28 | +4.8 % |
| Sn | 7.21 | −1.7 % | 7.38 | +0.6 % | 7.96 | +8.4 % |

- **Con −ε del ³P el error ya es menor que 2 % en Si, Ge y Sn.** La corrección CI (~0.2 eV) empuja
  todo hacia arriba y sobrestima, salvo en el estaño.
- **V_pol con el α_d calibrado empeora la ionización:** sobreliga al electrón de valencia. Es
  consistente con que α_d absorbe algo más que polarización (`PLAN-REGENERACION.md` §6).

### 2.6 Factor g del ³P₂ y mezcla con ¹D₂

| | g (HF+CI) | peso de ¹D₂ | g (CI+V_pol) | peso de ¹D₂ |
|---|---|---|---|---|
| Ge | 1.496164 | 0.77 % | 1.495316 | 0.94 % |
| Sn | 1.472668 | 5.47 % | 1.466527 | 6.69 % |

- **El 1.496 que la tesis atribuye al NIST** (hay que verificarlo) implica, a tres cifras, entre
  0.70 y 0.90 % de ¹D₂. Los dos cálculos del germanio caen dentro. La brecha de `discusion.tex:69`
  y la historia del CI que la cierra (`:71-73`) desaparecen.
- **Grado de mezcla:** el germanio no está «altamente mezclado» (menos de 1 %); el estaño sí
  mezcla de forma apreciable (5–7 %).

### 2.7 Otros números que cambian respecto a lo impreso

- **`tab:momentos_inversos`:** Ge da 4.793275; la tesis imprime 4.793285.
- **`tab:resultados_energia_global`:** el virial del Sn da 1.99999938 (la tesis imprime
  2.00000000), y la T del Ge, 2075.36257 (la tesis, 2075.36256).
- **`tab:integrales_slater`:**
  - Si: 0.32969971 y 0.16594819.
  - Ge: F⁰ = 0.31650266; su F² queda sin referencia.
  - Sn: con las cifras completas.
- **Pozo de V_rad del estaño:** −399.43 Ha en r = 0.042 a₀ (faltaba en `resultados.tex:98`).

---

## 3. Qué cambia en el manuscrito

Afirmaciones que los números contradicen, además de las de `FALTANTES-2026-09-13.md` §3.2:

| pasaje | dice hoy | dicen los números |
|---|---|---|
| `resultados.tex:248, 263-265`; `correlacion-estructura-fina.tex:163-167` | HF sobrestima el desdoblamiento y V_pol lo baja relajando el orbital hacia afuera | HF+CI lo subestima (−11 % Ge, −15 % Sn); V_pol contrae el orbital y sube ζ |
| `resultados.tex:292` | α_d = 3 acierta el ³P₂ y α_d = 4 el ³P₁; compromiso en 3.5 | los cruces son 2.87 y 3.22; α_d = 3.0 deja los dos a menos de 2 % |
| `resultados.tex:269-271` | V_pol apenas altera los singletes; el ¹S₀ del Ge conserva más de 3000 cm⁻¹ de error, que se cerraría ampliando el espacio CI | ζ es 3–5 veces más sensible que F², pero los singletes se mueven cientos de cm⁻¹; con V_pol el ¹S₀ del Ge queda a +0.6 % y el que falla es el ¹D₂ (+15 %); el espacio ya está convergido |
| `resultados.tex:186-208` | Koopmans empeora con Z y el CI corrige centésimas de eV | con −ε del ³P el error es menor que 2 % en Si–Sn; el CI corrige ~0.2 eV y sobrestima |
| `resultados.tex:294-311` | V_pol reduce el error de ionización del Ge a la mitad | V_pol lo empeora: 1.4 % → 4.8 % |
| `resultados.tex:267` | V_pol no se aplica a los ligeros porque su core es rígido | cierto, y además en el carbono ζ sobra: V_pol, que lo sube, lo empeoraría |
| `discusion.tex:11` | el carbono obedece fielmente la regla de intervalos de Landé | el modelo sí (R = 1.98); el experimento no (R = 1.65), por el espín-espín |
| `discusion.tex:14` | la desviación del silicio viene del promedio de configuración | venía del signo de f₂ y de una transcripción |
| `discusion.tex:19, 64-73` | Ge «altamente mezclado»; g = 1.4975 con una brecha que cierra un CI de excitaciones np² → nd² | Ge tiene 0.8 % de ¹D₂; g = 1.4962 ya es compatible con 1.496; el CI es de pareja general, con parejas no equivalentes |

Contenido nuevo que ya tiene sus números:

- **Convergencia del CI:** `convergencia_ci.tex`, `convergencia_ci_costo.tex` y `convergencia_ci.pdf`.
- **C-DIIS:** `cdiis.tex` y `cdiis_convergencia.pdf`.
- **Diagnóstico de ζ:** `zeta.tex`, `zeta_razon.pdf` y `valores_texto.md`.
- **Barridos de V_pol:** `barrido_vpol_Ge.tex`, `barrido_vpol_Sn.tex` y `barrido_vpol.pdf`.
- **Factor g:** `lande.tex`.

---

## 4. Decisiones pendientes

1. **Carbono con m > 20.**
   - (a) Ampliar solo su curva a m = 24 y 28 y reportar m = 28 en su tabla. Es del orden de 10 min
     y ~1 GB de memoria, con un cambio chico en `config.jl` y en las etapas 2–4 (tamaños por
     elemento).
   - (b) Reportar m = 20 y la extrapolación.
   - Recomendación: (a). La afirmación de que el ¹S₀ alto es del modelo y no del tamaño del espacio
     necesita un carbono convergido.
2. **Figura de niveles.** Hoy muestra solo HF+CI y NIST; agregar la columna CI+V_pol de Ge y Sn es
   opcional y cuesta minutos.
3. **Literatura que el código no puede dar:**
   - F²(4p,4p) de Froese Fischer para el germanio y el virial de referencia del estaño;
   - g del ³P₂ del NIST para C, Si y Sn, y confirmar el 1.496 del germanio;
   - versión y fecha de consulta del NIST ASD, y el potencial de ionización del Sn (7.34 eV);
   - si el ζ(2p) = 31.946 cm⁻¹ de Froese Fischer es el de Blume-Watson;
   - polarizabilidades dipolares de los cores de Ge y Sn, para compararlas con α_d = 0.75 y 3.0;
   - el factor de escala de los F^k en Cowan (1981).

---

## 5. Orden de redacción sugerido

1. **`resultados.tex`:** las tablas por `\input` (correspondencia en `PLAN-REGENERACION.md` §5) y
   después las secciones en el orden del capítulo, con la tabla de la sección 3.
2. **Subsecciones nuevas:** convergencia del CI y C-DIIS en resultados; el diagnóstico de ζ al
   final de resultados o en la discusión.
3. **Capítulos de método y apéndice:**
   - `ciclo-scf.tex`: `sec:ci_method`, la formulación de C-DIIS y la calibración de V_pol;
   - `correlacion-estructura-fina.tex`: el mecanismo de V_pol;
   - el Apéndice F.
4. **`discusion.tex`,** completo.
5. **Cierre:** `conclusiones.tex` con trabajo futuro, la introducción y la pasada final de
   `FALTANTES-2026-09-13.md` §3.3.
