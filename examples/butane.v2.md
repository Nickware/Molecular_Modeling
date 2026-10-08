# Dinámica molecular de butano líquido: propuesta v2

Este protocolo propone una simulación clásica de muchas moléculas de n-butano con GROMACS y una topología OPLS-AA obtenida de una fuente de parametrización, por ejemplo [LigParGen](https://ligpargen.com/). Es independiente del cálculo DFT de una molécula aislada en `butane/butane.py`.

Los valores de temperatura, presión, duración, tamaño de caja y opciones de corte que aparecen aquí son una **configuración inicial**, no resultados validados. Hay que registrar las versiones del software y de la topología, confirmar la fase termodinámica, comprobar convergencia y contrastar con datos experimentales antes de presentar propiedades como resultados.

## 1. Definir el sistema y las condiciones

Para estudiar un líquido, primero elige una temperatura y presión a las que n-butano sea líquido. Como punto de partida de este protocolo puede investigarse `T = 250 K`, `P = 1 bar`; confirma que ese estado está en la región líquida con el diagrama de fases o la presión de vapor experimental que vayas a usar. No uses `298 K, 1 bar` como condición de líquido: a esas condiciones el butano estable es gas.

Anota antes de empezar:

- Temperatura y presión objetivo, con unidades.
- Fuente experimental y estado de fase para cada propiedad con la que compararás.
- Dos cantidades de moléculas para el estudio de tamaño finito; por ejemplo, `N = 256` y `N = 512`.
- Al menos tres semillas independientes por tamaño para estimar variabilidad entre réplicas.
- Criterios de equilibrio, duración de producción y plan de análisis, antes de observar los resultados.

GROMACS usa `nm` para longitud, `ps` para tiempo, `K` para temperatura, `bar` para presión y `kJ/mol` para energía en sus entradas MD.

## 2. Instalar GROMACS

Instala una versión estable de GROMACS en el host o en el contenedor usado para los ejemplos. Por ejemplo, en Ubuntu puedes usar el paquete de la distribución:

```bash
sudo apt update
sudo apt install -y gromacs
```

Comprueba que el ejecutable funciona y registra su versión:

```bash
gmx --version
```

Si tu sistema requiere paralelismo MPI, GPU o una versión concreta, sigue la guía oficial de instalación de [GROMACS](https://manual.gromacs.org/current/install-guide/index.html) y conserva las opciones de compilación en el registro del cálculo.

## 3. Obtener y conservar la parametrización OPLS-AA

En LigParGen, solicita **n-butano neutro** (`CCCC`) con el campo de fuerza **OPLS-AA** y formato de salida para GROMACS. Conserva el archivo descargado completo y todos sus archivos de topología y coordenadas. No reconstruyas cargas, tipos atómicos, reglas de combinación o términos enlazados copiando tablas generales del repositorio: la parametrización debe corresponder a esta molécula y al formato seleccionado.

Organiza los archivos, por ejemplo:

```text
butane-v2/
  parameters/
    ligpargen-original.zip
    ligpargen-output/
  runs/
    N256/rep01/
    N256/rep02/
    N256/rep03/
    N512/rep01/
    ...
  analysis/
```

Crea las carpetas de trabajo y, si LigParGen entrega un archivo ZIP, extrae una copia conservando el ZIP original:

```bash
mkdir -p butane-v2/parameters/ligpargen-output \
  butane-v2/runs/N256/rep01 butane-v2/runs/N256/rep02 \
  butane-v2/runs/N256/rep03 butane-v2/runs/N512/rep01 \
  butane-v2/runs/N512/rep02 butane-v2/runs/N512/rep03 \
  butane-v2/analysis
unzip butane-v2/parameters/ligpargen-original.zip \
  -d butane-v2/parameters/ligpargen-output
find butane-v2/parameters/ligpargen-output -maxdepth 2 -type f -print
```

Si el archivo se descarga con otro nombre o formato, conserva ese nombre original y usa la herramienta adecuada para inspeccionar/descomprimirlo. Los nombres y rutas de las coordenadas de los comandos posteriores son ilustrativos: ajústalos a los archivos que LigParGen realmente generó, sin renombrar ni modificar el original.

Registra en un `README` o manifiesto de la parametrización:

- Servicio o referencia de parametrización, fecha y versión disponible.
- Identidad química/SMILES, carga total, campo de fuerza y formato exportado.
- Nombres de cada archivo original, procedencia y, si es posible, sus sumas SHA-256 (`sha256sum archivo`).
- Cualquier conversión hecha desde el paquete original.

Inspecciona la topología antes de usarla. Conserva las definiciones de tipos y los parámetros de enlaces, ángulos, torsiones, cargas, interacciones de Lennard-Jones, interacciones 1-4 y sus factores de escala tal como fueron generados. Comprueba que la coordenada de una molécula y su topología describen los mismos 14 átomos y el mismo orden de átomos.

Al preparar la topología del sistema, parte de los archivos generados por LigParGen y conserva sus `#include`, `[ defaults ]`, `[ atomtypes ]`, `[ nonbond_params ]` y `[ moleculetype ]` según estén organizados. No añadas otra definición de OPLS ni dupliques `[ defaults ]`. En la sección `[ molecules ]`, el nombre debe coincidir exactamente con el nombre de `[ moleculetype ]`; el número se cambiará al total que realmente se haya insertado. `gmx grompp -pp processed.top` permite revisar el archivo que resulta tras procesar los includes.

No continúes si GROMACS informa que faltan tipos, parámetros, átomos o que el total de átomos de la topología no coincide con la estructura.

## 4. Preparar cajas periódicas independientes

Usa la coordenada monomolecular exportada por LigParGen como plantilla. Estima la longitud de caja inicial desde una densidad líquida de referencia aproximada; después la caja se relajará en NPT. Para `N` moléculas de masa molar $M$ y densidad inicial $\rho$, la longitud cúbica es

$$L = \left(\frac{N M}{\rho N_A}\,10^{21}\right)^{1/3} \quad \text{nm},$$

con $M$ en g/mol, $\rho$ en g/cm$^3$ y $N_A$ en mol$^{-1}$. El factor convierte cm$^3$ a nm$^3$. Por ejemplo, `rho_guess = 0.60 g/cm^3` es solo una estimación para empaquetar; no se debe reportar como densidad simulada.

Este comando calcula un valor inicial de caja para 256 moléculas; cambia `N` y `rho_guess` para cada tamaño:

```bash
python3 - <<'PY'
import math

n_molecules = 256
molar_mass = 58.1222  # g/mol, n-butano
rho_guess = 0.60      # g/cm^3, solo estimación inicial
avogadro = 6.02214076e23
box_nm = (n_molecules * molar_mass * 1e21 / (rho_guess * avogadro)) ** (1 / 3)
print(f"Longitud inicial: {box_nm:.3f} nm")
PY
```

Desde un directorio de réplica, centra la molécula plantilla y crea una caja de 256 moléculas. Sustituye la longitud por el valor calculado y usa una semilla distinta en cada réplica:

```bash
gmx editconf -f ../parameters/ligpargen-output/butane.gro \
  -o butane-centered.gro -c

gmx insert-molecules -ci butane-centered.gro -nmol 256 \
  -box 3.45 3.45 3.45 -try 1000 -seed 41001 \
  -o box.gro
```

`gmx insert-molecules` coloca moléculas con orientaciones y posiciones aleatorias evitando los solapamientos definidos por sus radios. Comprueba el mensaje final: deben haberse insertado todas las moléculas solicitadas. Si no, aumenta `-try`, revisa la caja/plantilla y no pongas en `[ molecules ]` un número mayor que el insertado. Confirma además que `box.gro` contiene `14 * N` átomos y una celda cúbica con el lado esperado.

Repite el procedimiento para `N = 512` y para cada réplica, con directorios y semillas separados. La topología de cada réplica debe tener en `[ molecules ]` exactamente el número de moléculas presentes en su `box.gro`. Mantén fijos los parámetros físicos entre tamaños y réplicas.

## 5. Minimizar la energía

Crea `em.mdp` en cada directorio. Los parámetros de interacción deben ser consistentes con la topología LigParGen; los cortes siguientes son un punto de partida habitual, pero se deben contrastar con la parametrización exportada y la literatura del campo de fuerza:

```ini
integrator              = steep
emtol                   = 1000.0
emstep                  = 0.01
nsteps                  = 50000
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.2
coulombtype             = PME
rcoulomb                = 1.2
vdwtype                 = Cut-off
rvdw                    = 1.2
DispCorr                = EnerPres
pbc                     = xyz
constraints              = h-bonds
constraint_algorithm    = lincs
```

Antes de usar estos cortes, verifica que cada dimensión de la caja sea mayor que dos veces el mayor corte no enlazado. Si LigParGen recomienda otro tratamiento electrostático o de dispersión para su topología, usa ese protocolo documentado y aplícalo igual en todo el estudio.

Construye y ejecuta la minimización:

```bash
gmx grompp -f em.mdp -c box.gro -p topol.top -o em.tpr -pp processed.top
gmx mdrun -deffnm em
```

Comprueba el potencial final, las fuerzas y que no queden contactos atómicos extremos. Una minimización elimina contactos desfavorables; no demuestra que el líquido esté equilibrado.

## 6. Equilibrar temperatura y presión (NPT)

Usa primero un NPT de equilibración con `v-rescale` para temperatura y acoplamiento de presión Berendsen solo como etapa de relajación. Berendsen no produce las fluctuaciones correctas de un ensamble NPT y no se debe usar para promediar propiedades finales.

Ejemplo de `npt-eq.mdp` para 250 K y 1 bar; ajusta `ref_t`, `ref_p` y `gen_seed` para el estudio:

```ini
integrator              = md
dt                      = 0.002
nsteps                  = 250000
continuation            = no
constraints              = h-bonds
constraint_algorithm    = lincs
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.2
coulombtype             = PME
rcoulomb                = 1.2
vdwtype                 = Cut-off
rvdw                    = 1.2
DispCorr                = EnerPres
pbc                     = xyz
tcoupl                  = v-rescale
tc-grps                 = System
tau_t                   = 1.0
ref_t                   = 250
pcoupl                  = Berendsen
pcoupltype              = isotropic
tau_p                   = 5.0
ref_p                   = 1.0
compressibility         = 1.0e-4
gen_vel                 = yes
gen_temp                = 250
gen_seed                = 41001
nstxout-compressed      = 1000
nstenergy               = 1000
nstlog                  = 1000
```

`compressibility` es un valor inicial de líquido que debe revisarse con la fuente del modelo/propiedades experimentales; no es un parámetro OPLS universal. Corre la etapa:

```bash
gmx grompp -f npt-eq.mdp -c em.gro -r em.gro -p topol.top -o npt-eq.tpr
gmx mdrun -deffnm npt-eq
```

Examina temperatura, presión, densidad, volumen y energía. Extiende la equilibración si la densidad deriva o si los observables aún no son estacionarios. No decidas la duración solo por haber alcanzado el número de pasos del ejemplo.

### Producción NPT para densidad y energía

Después del equilibrio, corre producción NPT con un barostato que muestree el ensamble, por ejemplo Parrinello-Rahman, y conserva todos los parámetros de interacción. Guarda como `npt-prod.mdp`:

```ini
integrator              = md
dt                      = 0.002
nsteps                  = 2500000
continuation            = yes
constraints              = h-bonds
constraint_algorithm    = lincs
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.2
coulombtype             = PME
rcoulomb                = 1.2
vdwtype                 = Cut-off
rvdw                    = 1.2
DispCorr                = EnerPres
pbc                     = xyz
tcoupl                  = v-rescale
tc-grps                 = System
tau_t                   = 1.0
ref_t                   = 250
pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau_p                   = 5.0
ref_p                   = 1.0
compressibility         = 1.0e-4
gen_vel                 = no
nstxout-compressed      = 1000
nstenergy               = 1000
nstlog                  = 1000
```

Con `dt = 0.002 ps`, `nsteps = 2500000` corresponde a 5 ns. El ejemplo conserva el estado del checkpoint y no genera velocidades nuevas. Usa el mismo tratamiento de interacciones y las mismas condiciones objetivo que en la equilibración; extiende producción si los promedios por bloques no convergen. Mantén el registro de energía suficientemente frecuente para el análisis por bloques.

```bash
gmx grompp -f npt-prod.mdp -c npt-eq.gro -t npt-eq.cpt \
  -p topol.top -o npt-prod.tpr
gmx mdrun -deffnm npt-prod
```

El tiempo de producción es un punto de partida: debe cubrir varias escalas de autocorrelación de la densidad y la energía.

## 7. Producción NVT para difusión

Usa una configuración de la región equilibrada NPT como condición inicial; NVT mantiene fija la caja de ese snapshot. Guarda como `nvt-prod.mdp`:

```ini
integrator              = md
dt                      = 0.002
nsteps                  = 2500000
continuation            = yes
constraints              = h-bonds
constraint_algorithm    = lincs
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.2
coulombtype             = PME
rcoulomb                = 1.2
vdwtype                 = Cut-off
rvdw                    = 1.2
DispCorr                = EnerPres
pbc                     = xyz
tcoupl                  = Nose-Hoover
tc-grps                 = System
tau_t                   = 2.0
ref_t                   = 250
pcoupl                  = no
gen_vel                 = no
nstxout-compressed      = 1000
nstenergy               = 1000
nstlog                  = 1000
```

Con `dt = 0.002 ps`, la producción de ejemplo dura 5 ns. Usa la misma temperatura y los mismos potenciales/cortes que en NPT. Nose-Hoover evita el termostato estocástico durante la medida de una propiedad dinámica como difusión; registra el termostato y su constante de acoplamiento porque también forman parte del protocolo. Para comparar réplicas del mismo tamaño, selecciona snapshots de regiones equilibradas y registra sus volúmenes y densidades; no fuerces una caja común distinta del volumen equilibrado.

```bash
gmx grompp -f nvt-prod.mdp -c npt-prod.gro -t npt-prod.cpt \
  -p topol.top -o nvt-prod.tpr
gmx mdrun -deffnm nvt-prod
```

Guarda coordenadas sin solvente en intervalos suficientemente cortos para el MSD; por ejemplo, cada 1–2 ps al principio. Revisa el tamaño de trayectoria que esto produce y aumenta la duración hasta que haya una región lineal del MSD bien muestreada.

## 8. Analizar las propiedades

### Densidad y energía de NPT

Extrae densidad, potencial, energía total, temperatura y presión de la trayectoria NPT:

```bash
gmx energy -f npt-prod.edr -o npt-observables.xvg
```

GROMACS solicitará los términos que quieres exportar; selecciona `Density`, `Potential`, `Total-Energy`, `Temperature` y `Pressure` si están disponibles. Descarta el transitorio de equilibración. Reporta la media de producción y su incertidumbre, no el último valor instantáneo. Anota las unidades impresas por GROMACS.

### Difusión molecular desde el MSD

`gmx msd -mol` calcula el MSD de centros de masa moleculares y corrige moléculas partidas por las fronteras periódicas. Usa la producción NVT sin incluir equilibration:

```bash
gmx msd -f nvt-prod.xtc -s nvt-prod.tpr -mol \
  -o msd.xvg -tu ps -trestart 10
```

Selecciona el grupo que contiene todas las moléculas de butano. Inspecciona `msd.xvg` y ajusta el intervalo de regresión solo sobre la región lineal de tiempo largo, sin la región balística inicial ni la cola ruidosa con pocas referencias. En tres dimensiones, la relación de Einstein es

$$D = \frac{1}{6}\frac{d}{dt}\langle |\mathbf{r}(t)-\mathbf{r}(0)|^2\rangle.$$

Usa el coeficiente que reporta `gmx msd` y conserva su unidad tal como la imprime el programa. Reporta también el intervalo de ajuste. Si aplicas una corrección de tamaño finito como Yeh–Hummer, presenta tanto $D_L$ sin corregir como la corrección, indicando la viscosidad y sus unidades/fuente.

## 9. Réplicas, incertidumbre y tamaño finito

Para cada tamaño (`N = 256`, `512`, como punto de partida), genera al menos tres réplicas realmente independientes. Cada una requiere:

- Semilla única para la inserción y para las velocidades iniciales.
- Caja/configuración inicial independiente; cambiar solo el nombre de archivo no crea otra réplica.
- Equilibración y producción completas, con sus `.mdp`, `.tpr`, `.log`, `.edr`, `.xtc`, checkpoints y topología archivados.

Estima la autocorrelación de los observables y elige bloques más largos que su tiempo de correlación característico. Calcula medias por bloque, incertidumbre por réplica y dispersión entre réplicas. No trates cada frame de la trayectoria como una muestra independiente. Para una propiedad $X$, reporta para cada réplica media, error estándar por bloques y duración analizada; después reporta la media y dispersión entre réplicas. Compara estas estadísticas entre tamaños. La densidad debe converger dentro de las incertidumbres; una dependencia apreciable de $D$ con el lado de caja indica efecto de tamaño finito, no un resultado ya convergido.

## 10. Comparación y registro de resultados

Compara con mediciones experimentales de n-butano líquido a las mismas temperatura y presión. Registra cita/DOI, unidades, fase, conversión de unidades y cualquier corrección aplicada. No compares una simulación de líquido con el valor experimental de gas ni con datos a otra temperatura.

Entrega, como mínimo, una tabla con:

| N moléculas | Réplicas | T / P | Densidad media ± incertidumbre | Energía media ± incertidumbre | $D_L$ media ± incertidumbre | Datos experimentales y fuente |
|---:|---:|---:|---:|---:|---:|---|
| 256 | 3 o más | condiciones del estudio | por calcular | por calcular | por calcular | cita |
| 512 | 3 o más | mismas condiciones | por calcular | por calcular | por calcular | cita |

No rellenes la tabla con valores esperados como si fueran resultados: las cifras deben venir de las trayectorias analizadas.

Conserva versiones de GROMACS y del campo de fuerza, archivos originales de LigParGen, coordenadas, topologías procesadas, entradas `.mdp`, semillas, logs, scripts de análisis y datos experimentales citados. Así se puede reconstruir cada réplica sin depender de los valores predeterminados del equipo.
