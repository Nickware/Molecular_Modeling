# Cálculo de energía de una molécula de agua (H₂O)

Este script es un archivo de entrada para NWChem que realiza la optimización geométrica seguida del cálculo de energía de una molécula de agua (H₂O) utilizando el método DFT con el funcional híbrido B3LYP y la base 6-31G*, sobre una geometría inicial especificada en coordenadas cartesianas en angstroms.

### Descripción del flujo

- **Definición de la molécula:** La sección `geometry` especifica las coordenadas iniciales de los tres átomos (oxígeno y dos hidrógenos) en unidades de angstroms, con opciones para imprimir el archivo XYZ y aplicar simetría automática, facilitando la preparación de la estructura para el cálculo.
- **Configuración de memoria y basis set:** Se asigna 8000 MB de memoria y se selecciona la base 6-31G* de la biblioteca estándar de NWChem para todos los átomos, proporcionando un balance entre precisión y costo computacional adecuado para moléculas pequeñas.
- **Parámetros DFT:** La sección `dft` configura el funcional de intercambio-correlación híbrido B3LYP (multiplet 1 para estado singlete cerrado) y define la carga neutra del sistema (charge 0).
- **Tareas secuenciales:** Primero ejecuta `task dft optimize` para minimizar la energía geométrica ajustando las posiciones atómicas hasta convergencia, y luego `task dft energy` para calcular con precisión la energía del mínimo optimizado, permitiendo análisis termodinámico.

### Utilidad y comentarios

- Este input es ideal para cálculos rápidos y precisos de geometrías y energías de moléculas pequeñas en química cuántica, usando NWChem de forma directa sin interfaces como ASE, lo que ofrece control total sobre parámetros nativos del programa.
- La geometría inicial parece representar una configuración razonable de H₂O cerca de su mínimo, pero la optimización refinará enlaces O-H y ángulo H-O-H automáticamente con B3LYP/6-31G*, comúnmente usado para propiedades espectroscópicas y reactividad.
- Requiere NWChem instalado; es ejecutable directamente desde la línea de comandos (ej. `nwchem input.nw output.nw`), generando archivos de salida con energías, geometrías finales y gradientes.

### Recursos adicionales

- La documentación de NWChem detalla opciones avanzadas para `dft` (grid, convergencia), basis sets personalizados y análisis de frecuencias (agregando `task dft freq`).
- Manuales y ejemplos de NWChem permiten extender a solventes implícitos, spin no cerrado o métodos post-HF para estudios más complejos.

Este formato nativo de NWChem es altamente flexible para automatización en flujos de trabajo computacionales y validación contra métodos experimentales.
