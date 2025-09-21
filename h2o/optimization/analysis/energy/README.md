### Análisis y Visualización de Energía DFT durante la Optimización

Este script permite **extraer y graficar la convergencia de la energía total DFT** (Density Functional Theory) a partir del archivo de salida `optimizacion.out` generado por el cálculo computacional implementado en `Nwchem` . El objetivo es visualizar cómo evoluciona la energía a lo largo de los pasos de optimización geométrica.

#### Características

- **Extracción automática** de valores de energía total DFT desde el archivo de salida.
- **Conversión de unidades** de Hartree a electronvoltios (eV).
- **Visualización gráfica** de la convergencia energética durante la optimización.
- **Personalización** de la gráfica para facilitar su análisis.


### Uso

#### 1. Requisitos

- Python 3.x
- Bibliotecas: `matplotlib`

Instala la biblioteca necesaria si no la tienes:

```bash
pip install matplotlib
```


#### 2. Archivo de entrada

Asegúrarse de tener el archivo `optimizacion.out` en el mismo directorio. Este archivo debe contener líneas con el texto `"Total DFT energy"` seguido del valor de energía.

Ejemplo de línea válida en `optimizacion.out`:

```
Total DFT energy: -152.345678 Hartree
```


#### 3. Ejecución

Guarda el script como `plot_dft_energy.py` y ejecútalo en la terminal:

```bash
python plot_dft_energy.py
```

El script:

- Lee cada línea del archivo buscando `"Total DFT energy"`.
- Extrae el valor numérico de energía y lo almacena.
- Convierte las energías de Hartree a eV.
- Muestra una gráfica de la energía en función del paso de optimización.


### Ejemplo de Gráfica

La gráfica generada mostrará la **convergencia de la energía total DFT** en cada paso de optimización, permitiendo identificar cuándo el sistema alcanza su energía mínima.

### Notas

- Si el formato de las líneas de energía del archivo es diferente, ajustar la expresión regular en el script.
- Se puede modificar los parámetros de la gráfica para adaptarla a sus necesidades.


### Autor

Este script está diseñado para facilitar el análisis de resultados en optimización molecular usando DFT. Si tiene sugerencia o encuentra algún error, no dude en contribuir.

