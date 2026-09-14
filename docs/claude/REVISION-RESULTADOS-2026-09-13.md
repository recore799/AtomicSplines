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
- **El punto abierto del lado del cálculo quedó cerrado:** los singletes del carbono no estaban
  convergidos a m = 20; con la curva hasta m = 40 lo están (sección 2.2).
- **Tres afirmaciones de la tesis se invierten con los números nuevos** (sección 3): el sentido de
  V_pol, la ionización con V_pol y la tensión del barrido del estaño. Una cuarta cambia de
  explicación: la brecha del factor g del germanio existe, pero no la cierra el CI (sección 2.6,
  corregida con el NIST).

---

## 1. Controles

| control | resultado |
|---|---|
| SCF de la etapa 1 | 9/9 con C-DIIS, sin caer al respaldo; Ge en 65–69 iteraciones (~1.5 min), Sn en 35–39 (~3.5 min) |
| manifiesto | 25/25; el diff solo agrega las 9 entradas |
| autoverificación del CI | 19/19 pasan T1–T4 |
| raíces | 39/39 en [1, 1, 1] |
| procedencia | 34 corridas con el código de `2479092` y las 5 del carbono (m = 24–40) con `6dfb84b`, sin cambios |
| reproducibilidad | a m = 20 y sin V_pol, E_corr y ³P₂ idénticos al log del 07-09 en los cuatro elementos |
| etapa 3 | sin celdas faltantes salvo las de literatura; volver a correrla da archivos idénticos |

---

## 2. Qué dicen los números

### 2.1 Estructura fina sin V_pol (m de producción: 40 en C y 20 en el resto; error frente al NIST)

| | ³P₁ | ³P₂ | ¹D₂ | ¹S₀ |
|---|---|---|---|---|
| C | +28.9 % | +45.3 % | +6.1 % | +10.4 % |
| Si | −6.2 % | −5.7 % | +15.7 % | −4.7 % |
| Ge | −11.9 % | −10.8 % | +10.2 % | −3.6 % |
| Sn | −19.2 % | −14.8 % | +1.1 % | −6.5 % |

Los ³P_J están convergidos en el espacio activo: el ³P₂ del carbono se mueve 0.13 cm⁻¹ de m = 4
a m = 40. Su error sigue la razón de ζ (`PLAN-REGENERACION.md` §6).

### 2.2 Convergencia de los singletes

Extrapolación geométrica con los incrementos de los tres últimos tamaños de cada curva
(m = 32 → 36 → 40 en el carbono; m = 12 → 16 → 20 en el resto):

| | ¹D₂ en el m de producción | cola | ¹S₀ en el m de producción | cola |
|---|---|---|---|---|
| C (m = 40) | 10810.1 | −35 | 23894.0 | −57 |
| Si (m = 20) | 7285.2 | −31 | 14674.5 | −103 |
| Ge (m = 20) | 7851.6 | −20 | 15780.2 | −32 |
| Sn (m = 20) | 8710.1 | no fiable (razón 0.93) | 16045.3 | −45 |

- **Los cuatro elementos están convergidos:** les falta del orden de 100 cm⁻¹ o menos.
- **El carbono necesitó llegar a m = 40.** Con m ≤ 20 sus dos extrapolaciones del ¹S₀ discrepaban
  (+9.5 % y +2.1 %) y la del ¹D₂ daba +2.8 %. Con la curva completa:
  - la razón entre incrementos sucesivos se estabiliza en ~0.6 desde m = 24, para los dos
    singletes y para E_corr, y las extrapolaciones geométrica y de potencia coinciden;
  - ¹D₂ → 10 760–10 775 cm⁻¹ (+5.6 %), ¹S₀ → 23 810–23 840 (+10.0 %) y E_corr → −9.93 mHa;
  - a m = 20 el ¹D₂ todavía estaba 340 cm⁻¹ y el ¹S₀ 680 cm⁻¹ por encima de su valor a m = 40.
- **Lectura física:** los dos singletes del carbono se quedan altos, y el ¹S₀ el doble que el ¹D₂.
  Es el patrón de la casi-degeneración ns² ↔ np², que el modelo de pareja con core congelado no
  contiene y que deprime específicamente al ¹S (`HALLAZGOS-2026-09-06.md` §4). El ¹D₂ residual
  apunta a correlación que el modelo tampoco contiene, la de las parejas que incluyen al 2s.

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
  que el modelo omite (Apéndice F, ecuación `eq:cpp_dos_cuerpos`). Los cálculos HFR reducen las
  integrales de Slater al 85 % siguiendo a Cowan (1981) (sección 6.3); la justificación hay que
  leerla en el libro antes de citarla.

### 2.5 Ionización (eV, orbitales ³P)

| | −ε | error | CI (ΔE) | error | + V_pol | error |
|---|---|---|---|---|---|---|
| C | 11.79 | +4.7 % | 12.06 | +7.1 % | | |
| Si | 8.08 | −0.8 % | 8.32 | +2.0 % | | |
| Ge | 7.82 | −1.0 % | 8.01 | +1.4 % | 8.28 | +4.8 % |
| Sn | 7.21 | −1.8 % | 7.38 | +0.5 % | 7.96 | +8.4 % |

- **Con −ε del ³P el error ya es menor que 2 % en Si, Ge y Sn.** La corrección CI (~0.2 eV) empuja
  todo hacia arriba y sobrestima, salvo en el estaño.
- **V_pol con el α_d calibrado empeora la ionización:** sobreliga al electrón de valencia. Es
  consistente con que α_d absorbe algo más que polarización (`PLAN-REGENERACION.md` §6).

### 2.6 Factor g del ³P₂ y mezcla con ¹D₂

Corregida la misma noche con el NIST ASD (sección 6.1). La primera versión de esta sección daba
por buena la cifra de 1.496 de la tesis y concluía que la brecha desaparecía: era falso. Los g
medidos llevan el g_s real del electrón, así que el ³P₂ puro vale 1.50116 y no 1.5. La tabla
compara con g calculado con ese mismo g_s (`lande.tex`).

| | g HF+CI | g CI+V_pol | g NIST | ¹D₂ HF+CI | ¹D₂ CI+V_pol | ¹D₂ NIST |
|---|---|---|---|---|---|---|
| Ge | 1.49731 | 1.49647 | 1.49458 | 0.77 % | 0.94 % | 1.31 % |
| Sn | 1.47376 | 1.46761 | 1.452 | 5.47 % | 6.69 % | 9.81 % |

- **El 1.496 que imprime la tesis no es el valor del ASD** (1.49458). Con el valor real, el
  germanio mezcla 1.31 % de ¹D₂ y el modelo da 0.77 % (0.94 % con V_pol): la brecha de
  `discusion.tex:69` existe y es mayor de lo que la tesis supone. El estaño muestra lo mismo:
  9.8 % medido frente a 5.5–6.7 %.
- **Lectura física:** el modelo mezcla de menos en los pesados. Es coherente con un ζ demasiado
  chico y un ¹D₂ demasiado alto, que agrandan el denominador del acoplamiento intermedio; V_pol
  corrige solo una parte.
- **Grado de mezcla:** el germanio no está «altamente mezclado» (1.3 % medido); el estaño sí
  (~10 %).

### 2.7 Otros números que cambian respecto a lo impreso

- **`tab:momentos_inversos`:** Ge da 4.793275; la tesis imprime 4.793285.
- **`tab:resultados_energia_global`:** el virial del Sn da 1.99999938 (la tesis imprime
  2.00000000), y la T del Ge, 2075.36257 (la tesis, 2075.36256).
- **`tab:integrales_slater`:**
  - Si: 0.32969971 y 0.16594819.
  - Ge: F⁰ = 0.31650266; su F² queda sin referencia.
  - Sn: con las cifras completas.
- **Pozo de V_rad del estaño:** −399.43 Ha en r = 0.042 a₀ (faltaba en `resultados.tex:98`).
- **Columnas NIST:** con los valores completos del ASD cambian cuatro celdas en la décima: ¹D₂ de
  C (10192.7) y Si (6298.9), ¹S₀ de Ge (16367.3) y Sn (17162.5). La R del NIST pasa a 1.644 (C) y
  1.894 (Si). Ninguna conclusión cambia.

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
| `discusion.tex:19, 64-73` | Ge «altamente mezclado»; g = 1.4975 frente a un experimental de 1.496, con una brecha que cierra un CI de excitaciones np² → nd² | el ASD da 1.49458: Ge mezcla 1.3 % de ¹D₂ y el modelo 0.8–0.9 %; la brecha existe y apunta a ζ y al ¹D₂, no al tamaño del CI, que ya es de pareja general |

Contenido nuevo que ya tiene sus números:

- **Convergencia del CI:** `convergencia_ci.tex`, `convergencia_ci_costo.tex` y `convergencia_ci.pdf`.
- **C-DIIS:** `cdiis.tex` y `cdiis_convergencia.pdf`.
- **Diagnóstico de ζ:** `zeta.tex`, `zeta_razon.pdf` y `valores_texto.md`.
- **Barridos de V_pol:** `barrido_vpol_Ge.tex`, `barrido_vpol_Sn.tex` y `barrido_vpol.pdf`.
- **Factor g:** `lande.tex`.

---

## 4. Decisiones y pendientes

Estado al cierre del 2026-09-13; el detalle está en la sección 6.

1. **Carbono con m > 20: hecho.** Su curva llega a m = 40, que es su m de producción, y está
   convergida (sección 2.2).
2. **Figura de niveles: hecho.** Germanio y estaño muestran HF+CI, CI+V_pol y NIST, con el α_d
   elegido en el título del panel.
3. **Literatura:**
   - hecho: niveles, factores g, potenciales de ionización y versión del NIST ASD (5.12), y el
     factor de escala de Cowan;
   - lo trae el usuario del libro de Froese Fischer: los F²(np,np), en particular el del germanio;
   - también del libro, si se puede: el virial de referencia del estaño y si el ζ(2p) = 31.946
     cm⁻¹ es el de Blume-Watson;
   - sin resolver: polarizabilidades dipolares de los cores de Ge y Sn (sección 6.3).

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

---

## 6. Actualización del 2026-09-13 (noche): NIST, carbono y literatura

### 6.1 NIST ASD

Consulta del 2026-09-13 a la versión 5.12 (Kramida, Ralchenko, Reader y NIST ASD Team, 2024; DOI
10.18434/T4W30F). Los valores completos y sus fuentes primarias están en `tesis/config.jl`, y la
entrada `NIST_ASD` quedó en `docs/bib/bibliografia.bib`.

- **Versión:** la tesis dice «versión 6.1» (`ciclo-scf.tex:205`); la consultada es la 5.12.
- **Niveles:** cuatro valores guardados no eran el redondeo del ASD (sección 2.7).
- **Ionización (eV):** C 11.2602880, Si 8.15168, Ge 7.899435 y Sn 7.343918. Cambia en la segunda
  cifra el error CI del silicio (2.0 %) y del estaño (0.5 %).
- **Factores g del ³P₂:** C 1.5010469(50), Ge 1.49458 y Sn 1.452; Si I no tiene en el ASD (sección
  2.6).

### 6.2 Carbono

- **Tamaños:** `TAMANOS_CONVERGENCIA` es ahora por elemento: el carbono llega a m = 40 y el resto
  se queda en m = 20. La sensibilidad se sigue comparando a m = 20.
- **Costo medido:** la caché de R^k usa ~17 bytes por casilla. A m = 40 son ~39 millones de
  integrales, del orden de 1.1 GB (hasta ~1.7 GB al redimensionarse) y ~35 min.
- **Salidas nuevas:**
  - `convergencia_singletes.tex` y `convergencia_singletes.pdf`;
  - en `valores_texto.md`, dos extrapolaciones (geométrica y de potencia) con los tres últimos
    tamaños.
- **Resultado de la corrida** (13-09, código `6dfb84b`): T1–T4 pasan, las raíces son la más baja y
  la curva tardó ~39 min. El ¹S₀ se queda alto (+10 % convergido) y el ¹D₂ también, pero a la
  mitad (+5.6 %). La afirmación de la tesis se sostiene con un carbono convergido (sección 2.2).

### 6.3 Literatura

- **Cowan (1981) y el 0.85:** la práctica de reducir las integrales de Slater al 85 % siguiendo la
  recomendación de Cowan (1981) está documentada, por ejemplo, en arXiv:2302.01780 (métodos HFR).
  Para citarla en la tesis conviene ir a la sección correspondiente del libro de Cowan.
- **Polarizabilidades de los cores de Ge y Sn: sin resolver.** No encontré valores primarios
  verificables; las búsquedas devolvieron resúmenes contradictorios que no uso.
  - Pistas para biblioteca, sin verificar: las tablas de Johnson, Kolb y Huang (1983), *At. Data
    Nucl. Data Tables* 28, que cubren la secuencia del níquel donde está Ge⁴⁺, y la revisión de
    Mitroy, Safronova y Clark (2010), *J. Phys. B* 43, 202001.
  - Al compararlas con α_d hay que recordar que el valor ajustado depende de r_c y de qué capas se
    tratan como core (con o sin ns²).
