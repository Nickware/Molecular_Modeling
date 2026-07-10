# Etano (C₂H₆): Perspectiva desde el Modelamiento Molecular

## Estructura y Conformación

El etano es la molécula **paradigmática** para estudiar la barrera rotacional en alcanos. A diferencia del butano, que presenta múltiples mínimos conformacionales, el etano tiene una simplicidad que lo convierte en **sistema de referencia** para parametrizar potenciales de torsión.

### Geometría de Equilibrio
- **Distancia C–C**: ~1.535 Å
- **Distancia C–H**: ~1.09 Å
- **Ángulo H–C–C**: ~111.2° (ligera compresión respecto al tetraédrico ideal de 109.5°)
- **Ángulo dihedral H–C–C–H**: 60° (conformación **escalonada**, *staggered*)

### Barrera Rotacional
La rotación alrededor del enlace C–C presenta:

| Conformación | Ángulo Dihedral | Energía Relativa |
|--------------|-----------------|------------------|
| **Escalonada** (*staggered*) | 60°, 180°, 300° | 0 (mínimo) |
| **Eclipsada** (*eclipsed*) | 0°, 120°, 240° | ~2.9 kcal/mol (~12.1 kJ/mol) |

Esta barrera de ~3 kcal/mol surge de la **repulsión torsional** entre orbitales de enlace C–H adyacentes (interacción de hiperconjugación vs. repulsión de Pauli), no de repulsión estérica pura como se pensaba clásicamente.

---

## Descripción Cuántica: Origen de la Barrera

### Cálculos *Ab Initio*
Los métodos de química cuántica revelan que la barrera rotacional del etano es un fenómeno **predominantemente cinético-electrónico**:

- **HF/STO-3G**: ~3.1 kcal/mol (sobrestima ligeramente)
- **MP2/cc-pVTZ**: ~2.9 kcal/mol
- **CCSD(T)/CBS**: ~2.88 kcal/mol (valor de referencia)
- **DFT (B3LYP/6-31G\*\*)**: ~2.8 kcal/mol

El análisis de descomposición de energía (EDA, Energy Decomposition Analysis) muestra:
- **Repulsión de Pauli** (intercambio): dominante en la conformación eclipsada
- **Hiperconjugación** (donación σ→σ\*): favorece la conformación escalonada (~4–5 kcal/mol de estabilización)
- **Interacción electrostática**: contribución menor (~0.5 kcal/mol)

> **Insight clave**: la barrera no se debe a repulsión estérica H···H (los H están a 2.5 Å en la eclipsada, lejos del contacto de van der Waals ~2.4 Å), sino a la **maximización del solapamiento orbital** en la conformación escalonada.

---

## Campo de Fuerzas para Etano

### Parametrización del Potencial de Torsión

El potencial dihedral se expresa típicamente como serie de Fourier:

$$V(\phi) = \frac{V_3}{2}[1 + \cos(3\phi)] + \frac{V_6}{2}[1 + \cos(6\phi)] + ...$$

Para el etano, el término dominante es $V_3$:

| Force Field | $V_3$ (kcal/mol) | $V_6$ (kcal/mol) | Origen de datos |
|-------------|------------------|------------------|-----------------|
| **OPLS-AA** | 2.9 | 0.0 | Ajuste a propiedades termodinámicas líquidas |
| **CHARMM** | 2.9 | 0.0 | Espectroscopía IR/Raman |
| **AMBER (GAFF)** | 2.9 | 0.0 | Cálculos HF/6-31G\* |
| **MM3** | 2.8 | 0.0 | Ajuste global a hidrocarburos |
| **TraPPE-UA** | 3.0 | 0.0 | Propiedades de fase (líquido-vapor) |

La periodicidad 3-fold refleja la simetría $C_3$ de los grupos metilo: tres mínimos equivalentes por rotación completa.

### Modelos de Átomos Unidos (United-Atom, UA)

En simulaciones de fluidos, el etano se modela frecuentemente como **dos centros interactivos** (CH₃–CH₃), reduciendo el costo computacional:

- **TraPPE-UA**: $\sigma_{CH_3} = 3.75$ Å, $\epsilon_{CH_3} = 98.0$ K
- **Potencial de Buckingham**: $U(r) = A\exp(-Br) - C/r^6$ (usado en cristales moleculares)

---

## Propiedades Mecánicas desde el Modelamiento

### 1. **Fase Gaseosa**

#### Ecuación de Estado y Segundo Coeficiente Virial
$$B_2(T) = -2\pi N_A \int_0^\infty \left[e^{-\beta u(r)} - 1\right] r^2 dr$$

Para etano a 300 K: $B_2 \approx -180$ cm³/mol (experimental: -184 cm³/mol). Los modelos UA con LJ puro reproducen bien $B_2$; los modelos AA (all-atom) mejoran la precisión a bajas T.

#### Propiedades de Transporte
| Propiedad | Método | Valor (298 K, 1 atm) |
|-----------|--------|----------------------|
| Viscosidad $\eta$ | MD + Green-Kubo | ~92 μPa·s (exp: 93) |
| Difusividad $D$ | MD (gas diluido) | ~0.15 cm²/s |
| Conductividad térmica $\lambda$ | MD + fluctuaciones | ~18 mW/m·K |

### 2. **Fase Líquida**

El etano líquido (T < 305 K a 1 atm) es un fluido molecular simple:

| Propiedad | Simulación MD | Experimental |
|-----------|-------------|--------------|
| Densidad (184 K) | ~0.546 g/cm³ | 0.546 g/cm³ |
| Entalpía de vaporización | ~14.7 kJ/mol | 14.7 kJ/mol |
| Tensión superficial | ~16.3 mN/m | 16.3 mN/m |
| Compresibilidad isotérmica $\kappa_T$ | ~1.8 GPa⁻¹ | ~1.8 GPa⁻¹ |

#### Estructura Local: Funciones de Distribución Radial
$$g_{CC}(r), \quad g_{CH}(r), \quad g_{HH}(r)$$

- Primer pico de $g_{CC}(r)$: ~4.0 Å (contacto C···C de moléculas vecinas)
- Coordinación de primeros vecinos: ~12 (empaquetamiento cercano a esférico)
- Orientación preferencial: los ejes C–C tienden al **empaquetamiento paralelo** en el primer shell, aunque con desorden rotacional significativo.

### 3. **Fase Sólida**

El etano cristaliza en dos polimorfos:

| Fase | Estructura | Rango de T | Grupo Espacial |
|------|-----------|------------|----------------|
| **Fase I** | Cubica (plástica) | > 89.9 K | Fm3m |
| **Fase II** | Monoclínica | < 89.9 K | P2₁/c |

#### Fase Plástica (I)
- **Rotor libre**: las moléculas reorientan casi libremente mientras mantienen sus centros de masa en una red cúbica.
- **Entropía de transición**: ~6.4 J/mol·K (pequeña, consistente con orden orientacional parcial).
- Simulaciones MD con potenciales de orientación anisotrópicos (Gay-Berne, Kihara) reproducen la fase plástica.

#### Fase Monoclínica (II)
- **Orden orientacional completo**: las moléculas adoptan orientaciones definidas.
- Barrera de reorientación: ~2–3 kcal/mol (similar a la barrera torsional intramolecular, pero ahora intermolecular).
- Dinámica de salto de 60° entre pozos de potencial orientacional.

---

## Espectroscopía Vibracional y Dinámica

### Modos Normales (13 modos: 3N – 6 = 9 vibracionales + 3 rotacionales + 3 traslacionales)

| Modo | Frecuencia (cm⁻¹) | Actividad | Descripción |
|------|-------------------|-----------|-------------|
| $\nu_1$ (a₁g) | 2954 | Raman | Estiramiento C–H simétrico |
| $\nu_2$ (a₁g) | 1388 | Raman | Deformación simétrica (umbrella) |
| $\nu_3$ (a₁u) | 2896 | Inactivo | Estiramiento C–H simétrico (torsionado) |
| $\nu_4$ (e_g) | 2969 | Raman | Estiramiento C–H degenerado |
| $\nu_5$ (e_g) | 1379 | Raman | Deformación degenerada (rocking) |
| $\nu_6$ (e_u) | 2985 | IR | Estiramiento C–H degenerado |
| $\nu_7$ (e_u) | 1468 | IR | Deformación degenerada (scissoring) |
| $\nu_8$ (e_u) | 822 | IR | Deformación degenerada (wagging) |
| $\nu_9$ (a₂u) | 289 | IR | Torsión alrededor de C–C |

> La frecuencia de torsión $\nu_9$ (~289 cm⁻¹) está directamente relacionada con la barrera rotacional: $\nu_{tors} \propto \sqrt{V_3/I_{red}}$, donde $I_{red}$ es el momento de inercia reducido.

### Ensanchamiento de Líneas en Fase Condensada
En líquido/sólido, las colisiones inducen:
- **Dephasing vibracional** (T₂): ~1–2 ps
- **Reorientación molecular** ($\tau_{reor}$): ~2–5 ps en líquido a T cercana al punto de ebullición
- **Tiempo de correlación de torsión**: ~0.1 ps (libración en el pozo escalonado)

---

## Mecánica Estadística del Etano

### Función de Partición Configuracional
Para N moléculas en fase condensada:

$$Q = \frac{1}{N!h^{3N}} \int d\mathbf{p}^N d\mathbf{q}^N \exp\left[-\beta H(\mathbf{p},\mathbf{q})\right]$$

En la aproximación de separabilidad:
$$Q \approx q_{trans}^{3N} \cdot q_{rot}^{2N} \cdot q_{vib}^{9N} \cdot q_{tors}^{N} \cdot q_{inter}^{N(N-1)/2}$$

La contribución torsional:
$$q_{tors} = \int_0^{2\pi} e^{-\beta V(\phi)} d\phi \approx \sqrt{\frac{2\pi k_B T}{V_3}} \quad \text{(para } k_BT \ll V_3\text{)}$$

A temperatura ambiente ($k_BT \approx 0.6$ kcal/mol << $V_3$), el etano está **localizado** en los pozos escalonados, con saltos termicamente activados entre ellos.

### Transición de Fase Líquido-Vapor
- **Temperatura crítica**: $T_c = 305.32$ K, $P_c = 48.72$ bar, $\rho_c = 0.207$ g/cm³
- **Punto triple**: 89.89 K, 1.1 × 10⁻³ bar
- **Calor latente de vaporización**: $\Delta H_{vap} = 14.7$ kJ/mol (en el punto de ebullición normal, 184.6 K)

Las simulaciones de Monte Carlo en ensamble Gibbs (GEMC) o dinámica molecular en NPT con termostato de Nose-Hoover reproducen la curva de coexistencia con precisión del 2–5%.

---

## Comparación con Butano: Lecciones del Modelamiento

| Aspecto | Etano | Butano |
|---------|-------|--------|
| **Barreras conformacionales** | Una sola barrera (3-fold) | Múltiples mínimos (anti, gauche) |
| **Complejidad del PES** | Simple, periódico | Rico, con interacciones 1,4 |
| **Propiedades termodinámicas** | Fluido casi ideal | Desviaciones por conformacionalismo |
| **Fase sólida** | Plástica + monoclínica | Múltiples polimorfos (ortorrómbicos) |
| **Aplicación como benchmark** | Parametrización de torsión | Parametrización de no-enlace 1,4 |

El etano es el **"hidrógeno molecular"** de la química orgánica computacional: su simplicidad permite aislar efectos y validar metodologías antes de escalar a sistemas más complejos.

---

## Aplicaciones Avanzadas del Modelo

### 1. **Mezclas y Solubilidad**
- Etano + metano (gas natural): simulaciones GEMC predicen diagramas de fase.
- Etano + agua: formación de **hidratos de gas** (estructura sI). Las simulaciones de MD con potenciales polarizables predicen las condiciones de estabilidad.

### 2. **Adsorción en Materiales Porosos**
- Zeolitas, MOFs, carbones activados: GCMC con etano como adsorbato modela la separación C₂H₆/CH₄.
- Isotermas de adsorción a 298 K: buena concordancia con experimentos para modelos UA.

### 3. **Combustión y Cinética**
- Mecanismos de oxidación: reacciones elementales C₂H₆ + O₂ → C₂H₅ + HO₂ (barrera ~50 kcal/mol).
- Dinámica de trayectorias *ab initio* (AIMD) para estudiar la disociación térmica C₂H₆ → 2 CH₃ (D(C–C) ≈ 85 kcal/mol).

---

## Resumen del Paradigma

| Escala Espacial | Método | Propiedades Accesibles |
|-----------------|--------|------------------------|
| Electrónica (Å) | CCSD(T), DFT | Barrera torsional, espectro vibracional, disociación |
| Molecular (nm) | MD/FF (OPLS, TraPPE) | EoS, transporte, estructura local |
| Mesoscópica (μm) | DPD, LB | Flujo en poros, mezclas complejas |
| Continuo (cm) | CFD + EoS | Diseño de procesos, separaciones |

El etano encapsula la **filosofía del modelamiento molecular**: partir de la descripción cuántica precisa de una molécula simple, parametrizar un campo de fuerza transferible, y escalar a propiedades termodinámicas y de transporte de sistemas macroscópicos.

# Optimización geométrica de la molécula de etano utilizando la interfaz de ASE con NWChem

Este script realiza una optimización geométrica de la molécula de etano utilizando la interfaz de ASE (Atomic Simulation Environment) con el motor de cálculos cuánticos NWChem, aplicando el funcional de densidad PBE y un algoritmo de optimización BFGS sobre una estructura inicial almacenada en 'ethane.xyz'.[2][7]

### Descripción del flujo

- **Lectura y definición de la molécula:** El objeto `Atoms` se construye a partir del archivo 'ethane.xyz', que contiene la geometría inicial de la molécula en formato XYZ.[2]
- **Asignación del calculador:** Al objeto molecular se le asigna el calculador NWChem, usando como funcional el 'PBE'. Es posible definir más parámetros (basis set, convergencia, etc.), pero aquí se usa una configuración mínima. Los parámetros adicionales pueden especificarse como diccionarios en el objeto NWChem para personalizar aún más el cálculo.[7][2]
- **Optimización geométrica:** La clase BFGS realiza la minimización de la energía variando la geometría molecular hasta que la fuerza máxima sobre cualquier átomo sea menor que $0.01$ $eV/Å$ (criterio de convergencia).[2]
- **Almacenamiento y análisis:** Cuando la optimización concluye, la geometría final se guarda en 'ethane_opt.xyz' y se obtiene la energía potencial del sistema optimizado mediante `get_potential_energy()`.[2]

### Utilidad y comentarios

- Este script es una plantilla básica para realizar optimizaciones geométricas con NWChem dentro de ASE, útil tanto en estudios de química computacional como en investigación de materiales. Es esencial que NWChem esté correctamente instalado y configurado en el entorno donde se ejecuta ASE.[8]
- Resulta conveniente ajustar parámetros como el basis set, el tipo de funcional, criterios de convergencia y paths de los archivos según los objetivos del estudio y la infraestructura disponible.[7][2]
- El valor fmax puede adaptarse para mayor precisión o rapidez en la optimización según el caso.[2]

### Recursos adicionales

- La documentación oficial de ASE para el módulo NWChem detalla los parámetros opcionales y métodos disponibles para personalizar aún más este tipo de cálculos.[7][2]
- La documentación de NWChem y ejemplos de input permiten comprender cómo ASE traduce las definiciones de objetos Python a archivos de entrada para NWChem.[4][8]

Si se requieren explicaciones sobre la personalización de parámetros, la integración con scripts más complejos o la visualización de trayectorias de optimización, revisar la documentación asociada a ASE Python.

### Documentación adicional

[1](https://nwchemgit.github.io/Geometry-Optimization.html)
[2](https://wiki.fysik.dtu.dk/ase/ase/calculators/nwchem.html)
[3](https://www.youtube.com/watch?v=HEpZQMFi1LE)
[4](https://nwchemgit.github.io/Getting-Started.html)
[5](https://winmostar.com/en/tutorials/QM_tutorial_1(Basic).pdf)
[6](https://docs.hpc.taltech.ee/chemistry/nwchem.html)
[7](https://wiki.fysik.dtu.dk/ase/_modules/ase/calculators/nwchem.html)
[8](https://docs.rcc.uchicago.edu/software/apps-and-envs/nwchem/)
[9](https://www.ehu.eus/sgi/ARCHIVOS/Nwchem)
[10](http://www.fqt.izt.uam.mx/software_fqt/user_4.7/user/node38.html)
[11](https://www.youtube.com/playlist?list=PLkbV4NC7mZIXvouZ4vhWeQIEsw4GJ3v6I)
[12](https://www.theochem.ru.nl/quantumchemistry/2020/cp1.html)
[13](https://comp.chem.umn.edu/nwchemrate/190807_NWChemrate_19_Manual.pdf)
