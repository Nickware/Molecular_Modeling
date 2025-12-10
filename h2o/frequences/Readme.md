# Frecuencias vibracionales de una molécula de agua (H₂O)

Este script es un archivo de entrada para NWChem que realiza la optimización geométrica, cálculo de energía de un solo punto y análisis de frecuencias vibracionales de una molécula de agua (H₂O) utilizando DFT con el funcional híbrido B3LYP y la base 6-31G*, sobre una geometría inicial en coordenadas cartesianas en angstroms.

### Descripción del flujo

- **Definición de la molécula:** La sección `geometry` especifica las coordenadas iniciales de los tres átomos (oxígeno y dos hidrógenos) en angstroms, con opciones para generar archivo XYZ y simetría automática, preparando la estructura para las tres tareas secuenciales.
- **Configuración de memoria y basis set:** Se reserva 8000 MB de memoria y se usa la base 6-31G* de la biblioteca estándar para todos los átomos, balanceando precisión y eficiencia en cálculos de moléculas pequeñas.
- **Parámetros DFT:** La sección `dft` define el funcional B3LYP con multiplicidad 1 (estado singlete cerrado) y carga neutra (charge 0), adecuado para sistemas cerrados como H₂O en su estado fundamental.
- **Tareas secuenciales:** Ejecuta `task dft optimize` para minimizar la energía geométrica, `task dft energy` para refinar la energía en la geometría optimizada, y `task dft freq` para calcular el espectro vibracional completo (frecuencias, modos normales e intensidades IR), validando la naturaleza de mínimo (sin frecuencias imaginarias).

### Utilidad y comentarios

- Este input completo es ideal para caracterización termodinámica y espectroscópica de moléculas pequeñas, proporcionando geometría optimizada, energía precisa y espectro vibracional para comparación experimental directa.
- La adición de frecuencias permite calcular entalpías, entropías y energías libres a diferentes temperaturas mediante aproximaciones armónicas, esencial en cinética química y termodinámica.
- Se ejecuta directamente con NWChem (`nwchem input.nw output.nw`), generando salida detallada con geometrías, energías, matrices de fuerzas y espectro vibracional listo para visualización.

### Recursos adicionales

- La documentación de NWChem detalla opciones para `freq` (análisis Raman, isotópicos), solventación implícita y basis sets más grandes para mayor precisión.
- Los manuales incluyen ejemplos de automatización de múltiples conformeros y de análisis de transiciones de estado mediante optimizaciones con restricciones.

Este script nativo ofrece control total para estudios vibracionales completos y es base para workflows más avanzados en química computacional.
