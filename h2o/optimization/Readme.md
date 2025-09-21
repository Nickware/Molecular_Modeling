## Optimización de la estructura del Agua

Esta implementación esta diseñada para ejecutarse en un software de química computacional (como NWChem, Gaussian u ORCA) y aplica la **Teoría del Funcional de la Densidad** (DFT) para estudiar una molécula formada por oxígeno e hidrógeno. A continuación se resumen los elementos principales:

- **Definición de memoria, carga y título** para el cálculo.
- **Geometría molecular** especificada en Ångströms, imprimiendo la estructura y detectando su simetría.
- **Asignación de base** (`6-31G*`) para todos los átomos.
- **Configuración del método DFT**, usando el funcional híbrido B3LYP y multiplicidad 1 (estado singlete).
- **Tarea principal:** optimización de la geometría para encontrar la estructura de menor energía.


## Perspectivas bajo el mismo enfoque DFT

A partir del mismo esquema de script puedes solicitar diversas tareas cambiando la instrucción final (`task`). Algunas opciones frecuentes incluyen:

- **Optimización de geometría**
    - Encontrar la estructura más estable.
    - Ejemplo: `task dft optimize`
- **Cálculo de energía de punto único**
    - Calcula la energía electrónica y propiedades relacionadas sin cambiar la geometría.
    - Ejemplo: `task dft energy`
- **Análisis de frecuencias vibracionales**
    - Determinar las frecuencias normales e identificar mínimos o puntos de transición. Permite simular espectros IR/Raman.
    - Ejemplo: `task dft freq`
- **Propiedades electrónicas**
    - Obtención de propiedades como momento dipolar, cargas atómicas, distribución de la densidad electrónica, energías de orbitales HOMO-LUMO, etc.
    - Ejemplo: `task dft property`
- **Cálculo de transiciones electrónicas (TD-DFT)**
    - Simular espectros UV-Visible y análisis de excitaciones.
    - Ejemplo: `task tddft energy`
- **Cálculo de perfiles de reacción**
    - Escanear coordenadas internas para estudiar rutas de reacción o barreras energéticas.
- **Propiedades magnéticas**
    - Determinación de desplazamientos químicos para espectroscopía RMN.
- **Análisis de poblaciones y cargas**
    - Métodos como Mulliken, Löwdin o NBO pueden usarse si el programa lo permite.


### Ejemplo de sintaxis para distintas tareas

```bash
task dft optimize      # Optimización geométrica
task dft energy        # Energía de un solo punto
task dft freq          # Frecuencias vibracionales
task dft property      # Propiedades electrónicas varias
task tddft energy      # Transiciones electrónicas (TD-DFT)
```

Para implementar estas tareas bajo la misma metodología, cambiar la línea de la tarea al final del script, aprovechando toda la configuración previa bajo el enfoque DFT.
