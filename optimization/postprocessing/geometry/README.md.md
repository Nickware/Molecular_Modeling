<img src="https://r2cdn.perplexity.ai/pplx-full-logo-primary-dark%402x.png" class="logo" width="120"/>

# README.md

## Procesamiento y Análisis de Geometrías Moleculares en MATLAB

Este script en MATLAB permite **procesar, analizar y visualizar geometrías moleculares** extraídas de archivos de salida de NWChem. Está diseñado para automatizar la extracción de coordenadas atómicas, el cálculo de distancias interatómicas, ángulos, momentos de inercia y la visualización 3D de la última geometría optimizada.

### Características Principales

- **Extracción automática** de bloques de geometría desde archivos de salida de NWChem.
- **Cálculo de distancias** entre todos los pares de átomos.
- **Cálculo de ángulos** entre tríos de átomos.
- **Cálculo del tensor de momento de inercia** para cada geometría.
- **Visualización 3D** de la última geometría procesada, con etiquetas de elementos.


## Uso

### 1. Requisitos

- MATLAB (cualquier versión reciente)
- Archivo de salida de NWChem que contenga bloques de geometría en angstroms.


### 2. Estructura del Archivo de Entrada

El archivo debe contener bloques similares a:

```
geometry units angstroms
H    0.0000    0.0000    0.0000
O    0.0000    0.0000    0.9572
H    0.9266    0.0000   -0.2396
end
```


### 3. Ejecución

Guarda el script como `geometry_postprocessing.m` y ejecuta en MATLAB:

```matlab
geo_data = geometry_postprocessing('nombre_del_archivo.out');
```

- El script procesará todas las geometrías encontradas en el archivo y mostrará una visualización de la última geometría.
- La variable `geo_data` contendrá una celda con estructuras para cada geometría, incluyendo:
    - `atomos`: lista de símbolos de los átomos.
    - `coordenadas`: matriz de coordenadas (en Å).
    - `distancias`: matriz de distancias interatómicas.
    - `angulos`: lista de ángulos calculados.
    - `momento_inercia`: tensor de momento de inercia (3x3).


## Descripción de Funciones

| Función | Descripción |
| :-- | :-- |
| `geometry_postprocessing` | Función principal: procesa el archivo, calcula propiedades y visualiza. |
| `procesar_geometrias` | Extrae los bloques de geometría del archivo de salida. |
| `calcular_distancias` | Calcula todas las distancias entre pares de átomos. |
| `calcular_angulos` | Calcula todos los ángulos entre tríos de átomos. |
| `calcular_momento_inercia` | Calcula el tensor de momento de inercia para la geometría. |
| `visualizar_geometria` | Genera una visualización 3D de la geometría molecular. |
| `obtener_masas` | Devuelve las masas atómicas para los elementos soportados. |

## Ejemplo de Uso

```matlab
geo_data = geometry_postprocessing('output_nwchem.out');
```

- Se mostrará una figura 3D con la geometría y etiquetas de los átomos.
- Puedes acceder a las propiedades calculadas de cada geometría usando, por ejemplo: `geo_data{1}.distancias`.


## Notas

- Si tu archivo contiene elementos no listados en la función `obtener_masas`, deberás añadirlos manualmente.
- El script asume que los bloques de geometría están bien formateados y separados por la palabra clave `end`.
- Puedes modificar o extender las funciones para calcular otras propiedades moleculares según tus necesidades.


## Autor

Este script está pensado para facilitar el análisis estructural de moléculas optimizadas con NWChem, acelerando tareas comunes en química computacional y modelado molecular. Si tienes sugerencias o encuentras algún error, no dudes en contribuir o abrir un issue.

