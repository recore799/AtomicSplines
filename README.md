# AtomicSplines.jl

Librería en Julia para cálculos de estructura atómica ab-initio utilizando una base de **B-splines**.

Actualmente, el código se centra en la resolución del problema de estados ligados (Bound States) y el procedimiento de Campo Autoconsistente (SCF) para obtener configuraciones electrónicas precisas y estudiar su comportamiento bajo condiciones de confinamiento esférico.

## Resultados actuales

### 1. Hidrógeno bajo Confinamiento Esférico
Estudio del átomo de hidrógeno sometido a un potencial de confinamiento esférico impenetrable. El código permite analizar la evolución del sistema variando el radio de la caja ($R_c$).

La siguiente figura muestra:
* La evolución de las funciones de onda radiales.
* El aumento de la energía del estado base al reducir el radio de confinamiento.
* La presión ejercida por el átomo sobre la "pared" esférica en función del radio.

![Hidrógeno Confinado](hidrogeno_confinado.png)

### 2. Helio: Ciclo SCF y Apantallamiento
Resolución del átomo de Helio mediante el método de Hartree-Fock. Se ilustra la convergencia del ciclo SCF y se visualiza claramente el **efecto de apantallamiento** electrónico, donde la densidad electrónica efectiva modifica el potencial percibido.

![Efecto de Apantallamiento en Helio](helio_apantallamiento.png)

## Cómo usar

1. **Instalación:**
   ```julia
   git clone [https://github.com/tu-usuario/AtomicSplines.jl](https://github.com/tu-usuario/AtomicSplines.jl)
   cd AtomicSplines.jl
   julia --project=. -e 'using Pkg; Pkg.instantiate()'

2. **Ejecutar un script:**
    ```julia
    julia --project=. examples/03_hydrogen_pressure.jl
    
Actualmente, el repositorio cuenta con los siguientes scripts de ejemplo en la carpeta `examples/`:

* `examples/01_toy_poisson.jl`: **Validación FEM.** Resuelve la ecuación de Poisson para una onda sinoidal en 1D para verificar la implementación de Elementos Finitos.
* `examples/02_hydrogen.jl`: **Átomo de Hidrógeno.** Validación del método B-splines resolviendo la ecuación de Schrödinger radial.
* `examples/03_hydrogen_pressure.jl`: **Hidrógeno Confinado.** Exploración del confinamiento esférico, generando datos de energía y presión vs radio.
* `examples/04_helium.jl`: **Átomo de Helio.** Implementación del procedimiento SCF para sistemas multielectrónicos y visualización del efecto de apantallamiento.
