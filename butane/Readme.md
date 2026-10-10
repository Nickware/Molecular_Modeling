# Butano: optimización de geometría y análisis vibracional

Este ejemplo estudia una molécula aislada de butano (C4H10) mediante teoría del funcional de la densidad (DFT). El programa de [butane.py](butane.py) usa ASE como interfaz de cálculo y NWChem para evaluar la energía electrónica y las fuerzas. A partir de una geometría XYZ, optimiza las posiciones atómicas y calcula las frecuencias vibracionales en la aproximación armónica.

El ejemplo no simula butano líquido ni calcula propiedades de transporte o termodinámicas. Para instalar ASE y los calculadores usados por los ejemplos del repositorio, consulta la [guía general de instalación](../docs/installation/README.md).

## Física del problema

El butano tiene una cadena de cuatro carbonos y puede cambiar de conformación al rotar alrededor del enlace central C2–C3. La conformación anti suele ser la más estable; las dos conformaciones gauche tienen energía mayor, y las conformaciones eclipsadas presentan barreras torsionales. Esta diferencia se origina principalmente en interacciones torsionales y repulsiones entre átomos.

En este cálculo, la energía electrónica $E(\mathbf{R})$ se aproxima con DFT para posiciones nucleares $\mathbf{R}$ dadas. Bajo la aproximación de Born–Oppenheimer, los núcleos se optimizan sobre la superficie de energía potencial electrónica. Las fuerzas son $\mathbf{F}_i=-\nabla_i E$ y el optimizador BFGS cambia las coordenadas para llegar a una geometría estacionaria.

Las vibraciones se describen cerca de esa geometría mediante la matriz Hessiana, $H_{ij}=\partial^2E/\partial R_i\partial R_j$. ASE la estima mediante desplazamientos finitos y calcula los modos normales a partir de la Hessiana ponderada por las masas. Una curvatura negativa da una frecuencia imaginaria; su presencia significativa puede indicar que la estructura es un punto de silla y no un mínimo. Como el butano es una molécula no lineal de 14 átomos, se esperan $3N-6=36$ modos vibracionales.

## Metodología implementada

1. ASE lee las coordenadas iniciales de `butane.xyz`.
2. Se asigna el calculador NWChem con el funcional de intercambio y correlación PBE (`xc='PBE'`). Los demás detalles del método dependen de los valores predeterminados del calculador y de la instalación; el script no especifica explícitamente un conjunto de bases.
3. El optimizador BFGS relaja la geometría hasta que la fuerza máxima cae por debajo de `0.01 eV/Å` (`fmax=0.01`).
4. Se guarda la geometría optimizada y se imprimen la energía potencial y la fuerza máxima residual.
5. ASE calcula y resume las frecuencias vibracionales. El script cuenta las frecuencias cuya parte imaginaria supera `1e-6 cm^-1`, las representa como barras negativas y guarda una gráfica diagnóstica.

La clasificación de frecuencias debe leerse como una comprobación numérica, no como una prueba automática de la conformación global más estable. Una geometría inicial diferente puede converger a otro mínimo local; además, modos muy pequeños pueden ser sensibles a los parámetros numéricos y a la convergencia del cálculo.

## Archivos y resultados

| Archivo o salida | Descripción |
|---|---|
| `butane.xyz` | Geometría molecular de partida: 4 carbonos y 10 hidrógenos. |
| `optimization.traj` | Trayectoria de la optimización BFGS. |
| `butane_opt.xyz` | Geometría al terminar la optimización. |
| `diagnostico_frecuencias.png` | Gráfica de frecuencias y conteo de modos imaginarios. |
| Salida de consola | Energía potencial en eV, fuerza máxima residual en eV/Å y diagnóstico de modos imaginarios. |

ASE también puede crear archivos de caché para el cálculo vibracional en el directorio de trabajo. Los valores numéricos no están fijados aquí: dependen del método, de los valores predeterminados de NWChem/ASE, de la geometría inicial y de la convergencia alcanzada.

## Alcance y limitaciones

- El modelo es una molécula aislada de butano en fase gaseosa, sin moléculas vecinas ni entorno solvente. No se define una celda periódica y no se muestrean configuraciones de un ensamble estadístico. La optimización y el cálculo vibracional evalúan estructuras estáticas cerca de un mínimo; no integran ecuaciones de movimiento ni generan una trayectoria temporal de dinámica molecular.
- La optimización encuentra un punto estacionario cercano a la geometría inicial; no barre el ángulo torsional C–C ni garantiza encontrar el mínimo conformacional global.
- PBE sin corrección de dispersión explícita puede no describir con la misma precisión todas las contribuciones débiles que afectan diferencias conformacionales.
- El análisis vibracional es armónico alrededor de la geometría optimizada; no incluye anarmonicidad ni efectos térmicos explícitos.
- El umbral de `1e-6 cm^-1` aplicado a la parte imaginaria es muy estricto. Modos casi nulos y ruido numérico merecen inspección antes de concluir que hay un punto de silla.

## Perspectivas de mejora

1. Especificar explícitamente el conjunto de bases, carga, multiplicidad, tolerancias electrónicas y parámetros del cálculo para mejorar la reproducibilidad.
2. Generar conformaciones anti y gauche, optimizarlas por separado y comparar energías relativas; así se estudiaría la superficie conformacional que el README describe como contexto.
3. Añadir dispersión (por ejemplo, una corrección compatible con el método elegido) y contrastar las diferencias energéticas con un método o conjunto de bases de referencia.
4. Definir un umbral físico configurable para frecuencias imaginarias, reportar sus valores y validar la convergencia de la Hessiana antes de clasificar el mínimo.
5. Separar los archivos de salida por ejecución, registrar versiones de ASE/NWChem y guardar un resumen estructurado para facilitar comparación y trazabilidad.

## Referencias de ejecución

La [guía de entorno de instalación y configuración](../docs/installation/readme.md) documenta los requisitos, la instalación de NWChem y Python, y la ejecución del ejemplo en un contenedor Ubuntu administrado con Distrobox.