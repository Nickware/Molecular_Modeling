 # Pentano (C₅H₁₂): Perspectiva desde el Modelamiento Molecular

## Estructura e Isomería

El pentano es el **alcano más simple que exhibe isomería estructural**, lo que lo convierte en un sistema clave para validar la **transferibilidad** de parámetros de campo de fuerza.

### Tres Isómeros Estructurales

| Isómero | Fórmula | Simetría | Punto de Ebullición | Punto de Fusión |
|---------|---------|----------|---------------------|-----------------|
| **n-pentano** | CH₃–CH₂–CH₂–CH₂–CH₃ | C₂v (conformación anti-anti) | 309.2 K | 143.4 K |
| **Isopentano** (2-metilbutano) | (CH₃)₂CH–CH₂–CH₃ | C₁ (conformador más estable) | 301.0 K | 113.3 K |
| **Neopentano** (2,2-dimetilpropano) | C(CH₃)₄ | Td | 282.6 K | 256.6 K |

> **Insight clave**: la diferencia en puntos de ebullición (~27 K entre isómeros) es puramente **entrópica-estérica**. El neopentano, con su geometría esférica compacta, tiene menor área superficial y menor capacidad de empaquetamiento eficiente, reduciendo las interacciones intermoleculares.

---

## n-Pentano: Conformacionalismo y PES

### Grados de Libertad Torsionales

El n-pentano tiene **dos ángulos dihedrales** dominantes (C2–C3 y C3–C4), generando un **paisaje energético 2D**:

| Conformación | $\phi_{2345}$ | $\phi_{1234}$ | Energía Relativa (kcal/mol) | Población (~300 K) |
|--------------|---------------|---------------|----------------------------|-------------------|
| **tt** (trans-trans) | 180° | 180° | 0.00 | ~55% |
| **tg⁺** (trans-gauche) | 180° | +60° | ~0.5 | ~18% |
| **tg⁻** (trans-gauche) | 180° | –60° | ~0.5 | ~18% |
| **g⁺g⁻** (gauche-opuesto) | +60° | –60° | ~1.0 | ~4% |
| **g⁺g⁺** (gauche-sinclinal) | +60° | +60° | ~3.5 | ~0.1% |
| **g⁻g⁻** (gauche-sinclinal) | –60° | –60° | ~3.5 | ~0.1% |

### Interacción de Gauche-Pentano (Gauche-Pentane Interaction)

La conformación **g⁺g⁺** (y g⁻g⁻) presenta una **penalización severa** (~3.5 kcal/mol) debido a la **interacción de van der Waals 1,5** entre los extremos metilo:

$$\Delta E_{g^+g^+} \approx 2 \times E_{gauche} + E_{1,5-vdW}$$

Donde $E_{1,5-vdW}$ es la repulsión entre los carbonos C1 y C5, que en g⁺g⁺ se acercan a ~2.8 Å (menor que la suma de radios de van der Waals, ~3.4 Å).

### Potencial de Torsión en Campos de Fuerza

$$V(\phi_1, \phi_2) = V(\phi_1) + V(\phi_2) + V_{cross}(\phi_1, \phi_2)$$

El término cruzado $V_{cross}$ es **esencial** para reproducir la población de g⁺g⁺:

| Force Field | $V_{cross}$ incluido | Población g⁺g⁺ predicha | Precisión |
|-------------|---------------------|------------------------|-----------|
| OPLS-AA | No (aproximación aditiva) | ~2% | Subestima |
| CHARMM | No | ~1.5% | Subestima |
| AMBER (GAFF) | No | ~1.8% | Subestima |
| **TraPPE-UA** | Sí (implícito en ajuste global) | ~0.2% | Buena |
| **MM4** | Sí (término explícito) | ~0.1% | Excelente |

> **Lección del modelamiento**: la aproximación de **potenciales dihedrales independientes** falla para cadenas con más de 4 carbonos. Los términos cruzados (o ajustes globales) son necesarios para capturar la cooperatividad conformacional.

---

## Geometría de Equilibrio y Parámetros Estructurales

### n-Pentano (Conformación tt)

| Parámetro | Valor Experimental | B3LYP/6-31G** | MP2/cc-pVTZ |
|-----------|-------------------|---------------|-------------|
| r(C–C) | 1.534 Å | 1.531 Å | 1.533 Å |
| r(C–H) | 1.096 Å | 1.093 Å | 1.094 Å |
| ∠(C–C–C) | 113.5° | 113.2° | 113.4° |
| ∠(H–C–H) | 107.8° | 107.5° | 107.7° |
| ∠(H–C–C) | 111.0° | 110.8° | 110.9° |

La apertura del ángulo C–C–C respecto al tetraédrico (109.5°) compensa la **repulsión 1,3** entre hidrógenos en carbonos adyacentes.

### Isopentano (2-metilbutano)

El carbono ramificado (C2) presenta:
- Ángulo C–C–C central: ~111° (compresión por esterismo)
- Tres rotores metilo adicionales
- Conformación más estable: **anti** respecto al enlace C2–C3, con los metilos en posiciones escalonadas

### Neopentano

- Simetría Td perfecta (como metano pero con CH₃ en vez de H)
- Barrera de rotación del grupo metilo: ~4.5 kcal/mol (mayor que en etano por repulsión 1,5 H···H)
- Sin grados de libertad torsionales de cadena principal

---

## Campo de Fuerzas para Pentano

### Modelos de Átomos Unidos (UA)

| Modelo | Parámetros C–C | $\sigma_{CH_2}$ (Å) | $\epsilon_{CH_2}$ (K) | Ajuste |
|--------|---------------|---------------------|----------------------|--------|
| **TraPPE-UA** | Transferidos de butano | 3.93 | 46.0 | Propiedades de fase LV |
| **NERD** | Transferidos | 3.93 | 47.0 | Densidad, $\Delta H_{vap}$ |
| **SKS** | Optimizados para pentano | 3.905 | 47.5 | Segundo virial + líquido |

### Modelos All-Atom (AA)

| Modelo | Carga C | Carga H | $\sigma_{CC}$ (Å) | $\epsilon_{CC}$ (K) |
|--------|---------|---------|-------------------|---------------------|
| **OPLS-AA** | –0.12 | +0.06 | 3.50 | 33.0 |
| **CHARMM** | –0.18 | +0.06 | 3.40 | 30.0 |
| **AMBER (GAFF)** | –0.12 | +0.06 | 3.40 | 35.0 |

### Términos de No-Enlace Específicos para Pentano

Los campos de fuerza deben incluir:

1. **1,4-interacciones** (C1···C4, H···H): escala de 1/2 en OPLS-AA
2. **1,5-interacciones** (C1···C5): escala de 1/1 (completa) en algunos campos
3. **Términos de torsión cruzados**: $\phi_{2345} \times \phi_{1234}$

$$V_{cross} = K_{cross} \cdot [1 + \cos(3\phi_1)] \cdot [1 + \cos(3\phi_2)]$$

---

## Propiedades Mecánicas desde el Modelamiento

### 1. Fase Líquida

| Propiedad | Exp. (298 K, 1 atm) | MD TraPPE-UA | MD OPLS-AA |
|-----------|---------------------|--------------|------------|
| Densidad (g/cm³) | 0.626 | 0.626 | 0.631 |
| $\Delta H_{vap}$ (kJ/mol) | 26.4 | 26.2 | 26.8 |
| Tensión superficial (mN/m) | 15.5 | 15.2 | 15.8 |
| Viscosidad (mPa·s) | 0.22 | 0.19 | 0.24 |
| Difusividad (10⁻⁵ cm²/s) | 3.5 | 3.8 | 3.2 |
| Compresibilidad $\kappa_T$ (GPa⁻¹) | 2.15 | 2.10 | 2.20 |

#### Estructura Local

$$g_{CC}(r), \quad g_{CH}(r), \quad g_{HH}(r)$$

- **Primer pico $g_{CC}$**: ~5.2 Å (distancia centro-centro de moléculas vecinas)
- **Coordinación**: ~10–11 (menor que metano/etano por forma alargada)
- **Orientación preferencial**: en el primer shell, las cadenas tienden al **empaquetamiento paralelo** (nesting), maximizando el contacto vdW

#### Orden Orientacional

El parámetro de orden de Steinhardt–Nelson:

$$Q_l = \left(\frac{4\pi}{2l+1} \sum_{m=-l}^{l} |\bar{q}_{lm}|^2\right)^{1/2}$$

Para pentano líquido a 298 K:
- $Q_2 \approx 0.15$ (ligero orden nemático local)
- $Q_6 \approx 0.35$ (orden de empaquetamiento)

> A diferencia de alcanos más cortos, el pentano líquido desarrolla **dominios de orden parcial** (cybotactic clusters) donde las cadenas se alinean paralelamente sobre distancias de ~10–20 Å.

### 2. Fase Gaseosa

#### Segundo Coeficiente Virial
| T (K) | $B_2$ (cm³/mol) Exp. | TraPPE-UA | OPLS-AA |
|-------|---------------------|-----------|---------|
| 300 | –737 | –720 | –750 |
| 400 | –357 | –345 | –365 |
| 500 | –192 | –185 | –198 |

La forma alargada del pentano reduce $B_2$ (más negativo) respecto a metano/etano, indicando mayor **atracción efectiva** por mayor área superficial.

### 3. Fase Sólida

El n-pentano cristaliza en una fase **ortorrómbica** a baja temperatura:

| Parámetro | Valor |
|-----------|-------|
| Grupo espacial | P2₁2₁2₁ |
| Moléculas por celda | 4 |
| Conformación en sólido | tt (planar zig-zag) |
| Densidad (sólido, 100 K) | ~0.82 g/cm³ |
| $\Delta H_{fus}$ | 8.4 kJ/mol |

Las simulaciones de MD en NPT con celdas variables predicen correctamente la densidad del sólido, pero la **predicción de la estructura cristalina** requiere técnicas de predicción de estructura (USPEX, genetic algorithms, metadynamics).

---

## Espectroscopía Vibracional

### Modos Normales del n-Pentano (Conformación tt, C₂h)

Con 17 átomos: 3N – 6 = **45 modos vibracionales**:

| Región (cm⁻¹) | Modos | Actividad | Descripción |
|---------------|-------|-----------|-------------|
| 2800–3000 | 12 | Raman + IR | Stretching C–H simétrico/asimétrico |
| 1400–1500 | 8 | Raman + IR | Deformación C–H (scissoring, rocking) |
| 1300–1400 | 6 | Raman + IR | Wagging C–H |
| 1000–1200 | 5 | Raman + IR | Twisting C–H + stretching C–C |
| 800–1000 | 4 | Raman + IR | Rocking C–H + breathing |
| 200–500 | 6 | Raman | Torsión C–C (barreras conformacionales) |
| < 200 | 4 | Raman | Modos de red (libración, traslación) |

### Modos Torsionales Clave

| Modo | Frecuencia (cm⁻¹) | Descripción |
|------|-------------------|-------------|
| $\tau_1$ (C2–C3) | ~120 | Torsión interna, barrera ~3.5 kcal/mol |
| $\tau_2$ (C3–C4) | ~120 | Torsión interna, degenerado con $\tau_1$ |
| $\tau_{metilo}$ | ~250 | Rotación metilo terminal, barrera ~2.9 kcal/mol |

En fase condensada, los modos torsionales se **ensanchan** y **desplazan** debido a la inhibición de la rotación libre por el entorno.

---

## Comparación Conformacional: n-Pentano vs. Butano vs. Hexano

| Propiedad | Butano | Pentano | Hexano |
|-----------|--------|---------|--------|
| **Mínimos conformacionales** | 2 (anti, gauche) | 6 (tt, tg, gt, g⁺g⁻, g⁺g⁺, g⁻g⁻) | 9 (3 dihedrales) |
| **Barrera anti-gauche** | ~0.9 kcal/mol | ~0.5 kcal/mol | ~0.5 kcal/mol |
| **Penalización g⁺g⁺** | No aplica | ~3.5 kcal/mol | ~3.5 kcal/mol (por dihedral) |
| **Longitud de persistencia** (líquido) | ~4.5 Å | ~5.5 Å | ~6.5 Å |
| **Radio de giro** ($R_g$) | ~2.1 Å | ~2.6 Å | ~3.1 Å |

### Longitud de Persistencia

En el límite de cadena larga (polímero):

$$l_p = \frac{\langle R_{end}^2 \rangle}{2L_c}$$

Donde $L_c$ es la longitud de contorno. Para alcanos cortos:
- Pentano: $l_p \approx 5.5$ Å (comparable a su longitud total ~6.4 Å)
- Esto indica que el pentano es **semirrígido**: no completamente flexible, pero tampoco rígido.

---

## Isopentano y Neopentano: Efectos de Ramificación

### Isopentano (2-metilbutano)

#### Conformaciones
El grupo metilo en C2 introduce **asimetría**:

| Conformación | Energía Relativa | Población (~300 K) |
|--------------|-----------------|-------------------|
| Anti (metilos opuestos) | 0.00 | ~65% |
| Gauche (metilos adyacentes) | ~0.8 kcal/mol | ~35% |

La ramificación reduce la **longitud efectiva** de la molécula, disminuyendo la viscosidad y el punto de ebullición respecto al n-pentano.

### Neopentano (2,2-dimetilpropano)

#### Propiedades Únicas
- **Forma casi esférica**: radio de giro mínimo
- **Sin grados de libertad torsionales internos** (solo rotación de metilos)
- **Fase plástica**: rota libremente en fase sólida hasta muy baja T
- **Punto de fusión anómalo**: 256.6 K (mucho más alto que n-pentano: 143.4 K)

| Propiedad | n-pentano | Isopentano | Neopentano |
|-----------|-----------|------------|------------|
| $\Delta H_{vap}$ (kJ/mol) | 26.4 | 24.7 | 22.4 |
| $\Delta H_{fus}$ (kJ/mol) | 8.4 | 5.2 | 3.3 |
| $T_c$ (K) | 469.7 | 460.4 | 433.8 |
| $P_c$ (bar) | 33.7 | 33.8 | 32.0 |
| Factor acéntrico $\omega$ | 0.252 | 0.230 | 0.196 |

> **Lección de modelamiento**: la ramificación reduce el factor acéntrico y las interacciones intermoleculares efectivas, haciendo que el neopentano se comporte más como un **fluido casi esférico** (proto-Lennard-Jones).

---

## Mecánica Estadística y Función de Partición

### Contribución Conformacional

Para el n-pentano:

$$Q_{conf} = \int_0^{2\pi} \int_0^{2\pi} e^{-\beta V(\phi_1, \phi_2)} d\phi_1 d\phi_2$$

Aproximación discreta (suma sobre mínimos):

$$Q_{conf} \approx \sum_{i,j} g_{ij} \cdot e^{-\beta E_{ij}}$$

Donde $g_{ij}$ es la degenerancia (1 para tt, 2 para tg⁺/tg⁻, etc.).

### Entropía Configuracional

$$S_{conf} = R \ln Q_{conf} + RT \frac{\partial \ln Q_{conf}}{\partial T}$$

| Isómero | $S_{conf}$ (J/mol·K, 298 K) | Origen |
|---------|---------------------------|--------|
| n-pentano | ~18 | 6 conformaciones accesibles |
| Isopentano | ~15 | 2 conformaciones + 3 rotores metilo |
| Neopentano | ~8 | 4 rotores metilo (simetría restringe) |

---

## Simulaciones de Dinámica Molecular: Desafíos Específicos

### 1. Tiempo de Correlación Conformacional

$$\tau_{conf} = \int_0^\infty \frac{\langle \phi(0)\phi(t) \rangle}{\langle \phi^2 \rangle} dt$$

| Condición | $\tau_{conf}$ (ps) | Mecanismo |
|-----------|-------------------|-----------|
| Gas (1 atm, 400 K) | ~0.5 | Colisiones binarias |
| Líquido (298 K) | ~5–10 | Reorganización cooperativa de vecinos |
| Sólido (100 K) | > 1 ns | Saltos activados termicamente |

### 2. Difusividad y Mecanismo de Transporte

$$D = \frac{1}{6} \lim_{t\to\infty} \frac{d}{dt} \langle |\mathbf{r}(t) - \mathbf{r}(0)|^2 \rangle$$

| Fase | $D$ (10⁻⁵ cm²/s) | Mecanismo |
|------|-----------------|-----------|
| Gas (400 K, 1 atm) | ~50 | Movimiento balístico + colisiones |
| Líquido (298 K) | ~3.5 | Reptación + saltos |
| Sólido (100 K) | ~10⁻⁸ | Difusión intersticial (muy lenta) |

En líquido, el pentano exhibe un mecanismo de **reptación** (snake-like motion): la cadena se desliza a través de los vecinos manteniendo su conformación local.

### 3. Viscosidad: Efecto de la Longitud de Cadena

La viscosidad de alcanos líquidos sigue aproximadamente:

$$\eta \propto M^{3.4} \quad \text{(para polímeros en líquido)}$$

Para alcanos cortos:
- Metano: ~0.1 mPa·s
- Pentano: ~0.22 mPa·s
- Decano: ~0.86 mPa·s

Los modelos UA subestiman típicamente la viscosidad en ~10–20% porque no capturan la **fricción rotacional** completa de los grupos metileno.

---

## Comparación con la Serie de Alcanos Lineales

| Propiedad | Butano | Pentano | Hexano | Heptano |
|-----------|--------|---------|--------|---------|
| **N° de rotores** | 1 | 2 | 3 | 4 |
| **Complejidad PES** | 2D (1 dihedral) | 2D (2 dihedrales) | 3D (3 dihedrales) | 4D |
| **Población tt** | ~60% | ~55% | ~50% | ~45% |
| **Longitud de persistencia** | ~4.5 Å | ~5.5 Å | ~6.5 Å | ~7.5 Å |
| **Temperatura crítica** | 425 K | 470 K | 507 K | 540 K |
| **Modelo requerido** | Torsión simple | Torsión + cruzados | Torsión + cruzados + cooperatividad | Campo de fuerza de cadena |

---

## Métodos de Cálculo Avanzados

### Dinámica *Ab Initio* para Pentano

| Propiedad | Método | Resultado |
|-----------|--------|-----------|
| Barrera anti-gauche | MP2/cc-pVTZ | ~0.6 kcal/mol |
| Penalización g⁺g⁺ | CCSD(T)/CBS | ~3.4 kcal/mol |
| Barrera rotación metilo | B3LYP/6-311G** | ~2.9 kcal/mol |
| Frecuencias vibracionales | DFT (B3LYP) | Error < 2% respecto a exp. |

### Machine Learning para Alcanos

- **ANI-1x, ANI-2x**: potenciales neurales entrenados en DFT (ωB97X) que predicen energías de conformaciones de pentano con error ~0.5 kcal/mol
- **TorchMD**: simulaciones de MD con potenciales ML, velocidad ~100× más rápida que DFT
- **Aplicación**: predicción de propiedades de transporte de mezclas de alcanos

---

## Resumen del Paradigma de Modelamiento

| Escala | Método | Propiedades | Desafíos Específicos del Pentano |
|--------|--------|-------------|----------------------------------|
| Electrónica (Å) | DFT, MP2, CCSD(T) | Geometría, barreras, espectro | Cooperatividad de 2 dihedrales |
| Molecular (nm) | MD/FF (OPLS, TraPPE) | EoS, transporte, fases | Términos cruzados, 1,5-interacciones |
| Mesoscópica (μm) | DPD, MesoDyn | Mezclas, microemulsiones | Longitud de persistencia |
| Continuo (cm) | CFD + EoS (PC-SAFT) | Ingeniería de procesos | Isomería, ramificación |

El pentano representa la **frontera de la simplicidad** en el modelamiento molecular de alcanos:
- Con **2 rotores torsionales**, ya no es aditivo (como el butano)
- Con **3 isómeros**, prueba la transferibilidad de parámetros
- Con **forma semirrígida**, introduce efectos de orden orientacional en líquido

> **Insight final**: el pentano es el "puente" entre los alcanos cortos (metano, etano, butano), donde cada molécula es un caso especial, y los alcanos largos (decano, eicosano), donde las propiedades emergen de la estadística de cadenas flexibles. Su modelamiento exitoso requiere capturar tanto la **química conformacional** (barreras de torsión, interacciones 1,5) como la **física de fluidos** (empaquetamiento de cadenas, orden local).
