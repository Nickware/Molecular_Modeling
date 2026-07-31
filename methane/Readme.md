 # Metano (CH₄): Perspectiva desde el Modelamiento Molecular

## Estructura y Simetría

El metano es la **molécula de referencia** para la química orgánica computacional: la unidad estructural más simple con hibridación sp³ y simetría tetraédrica perfecta.

### Geometría de Equilibrio
| Parámetro | Valor Experimental | CCSD(T)/CBS |
|-----------|-------------------|-------------|
| Distancia C–H | 1.087 Å | 1.086 Å |
| Ángulo H–C–H | 109.47° | 109.47° |
| Grupo puntual | Td | Td |

### Propiedades Electrónicas
| Propiedad | Valor | Significado |
|-----------|-------|-------------|
| Momento dipolar $\mu$ | 0 D | Simetría tetraédrica anula momentos |
| Polarizabilidad $\alpha$ | 2.60 × 10⁻³⁰ m³ | Isotrópica (Td) |
| Cuadrupolo $\Theta$ | –2.6 × 10⁻²⁶ esu·cm² | Dominante en interacciones intermoleculares |
| Carga parcial C (NPA) | –0.54 e | Acumulación electrónica en el carbono |
| Carga parcial H (NPA) | +0.135 e | Polarización Cδ⁻–Hδ⁺ |

> **Insight clave**: aunque el metano es apolar (dipolo cero), su **cuadrupolo** no despreciable y su **polarizabilidad** dominan las interacciones intermoleculares. Los modelos que ignoran el cuadrupolo subestiman las propiedades de fase condensada.

---

## Descripción Cuántica: Hibridación y Orbitales

### Cálculos *Ab Initio*
| Método | Energía de Disociación C–H (kcal/mol) | Desviación geométrica |
|--------|--------------------------------------|----------------------|
| HF/6-31G** | 376 | ~0.01 Å |
| MP2/cc-pVTZ | 412 | < 0.005 Å |
| CCSD(T)/CBS | 419.3 | Referencia |
| DFT (B3LYP) | 405 | ~0.008 Å |
| DFT (PBE0) | 415 | ~0.006 Å |

### Análisis de Orbitales
- **Hibridación sp³ del carbono**: cuatro orbitales equivalentes dirigidos hacia los vértices de un tetraedro
- **Enlace C–H**: solapamiento sp³(C)–1s(H)
- **Orbitales moleculares**: 
  - 1a₁ (1s C, core)
  - 2a₁, 1t₂ (enlaces C–H)
  - 3a₁*, 2t₂* (antienlaces)
- **HOMO**: orbital triplemente degenerado 1t₂ (no enlazante, ligeramente C–H enlazante)
- **LUMO**: 3a₁* (antienlazante C–H)

La ausencia de pares solitarios en el carbono y la alta simetría hacen que el metano sea **químicamente inerte** y **espectroscópicamente simple**.

---

## Campo de Fuerzas para Metano

### 1. Modelo de Átomos Unidos (United-Atom, UA)

El metano es el sistema **ideal para UA**: un único centro de interacción con masa 16 u.

| Modelo | $\sigma$ (Å) | $\epsilon/k_B$ (K) | Origen |
|--------|-------------|-------------------|--------|
| **TraPPE-UA** | 3.730 | 148.0 | Propiedades de fase LV |
| **OPLS-UA** | 3.730 | 148.0 | Propiedades termodinámicas |
| **GROMOS** | 3.664 | 134.8 | Biomoléculas |
| **Mie (SAFT-γ Mie)** | 3.741 | 153.4 | EoS SAFT |

Potencial:
$$U(r) = 4\epsilon\left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right]$$

### 2. Modelos de Sitios Múltiples (All-Atom, AA)

Para mayor precisión en interacciones específicas:

| Modelo | Sitios | Características |
|--------|--------|----------------|
| **5-sites** | C + 4H + LJ en C | Cargas en C y H |
| **MP2-fit** | 5-sites | Parametrizado a dimeros MP2 |
| **Buckingham** | 5-sites | Repulsión exponencial |

Potencial 5-sites:
$$U = \sum_{i<j} \frac{q_i q_j}{r_{ij}} + 4\epsilon_{CC}\left[\left(\frac{\sigma_{CC}}{r_{CC}}\right)^{12} - \left(\frac{\sigma_{CC}}{r_{CC}}\right)^{6}\right] + \text{LJ cruzados}$$

### 3. Modelos con Cuadrupolo Explícito

Para capturar la anisotropía de la interacción:

$$U_{quad} = \frac{3}{4r^5}\left[\Theta_i : \Theta_j - 5(\Theta_i \cdot \hat{r})(\Theta_j \cdot \hat{r}) + \frac{35}{3}(\hat{r} \cdot \Theta_i \cdot \hat{r})(\hat{r} \cdot \Theta_j \cdot \hat{r})\right]$$

Los modelos con cuadrupolo mejoran la predicción del **segundo coeficiente virial** y la **tensión superficial**.

---

## Propiedades Mecánicas desde el Modelamiento

### 1. Fase Gaseosa

#### Ecuación de Estado y Segundo Coeficiente Virial
$$B_2(T) = -2\pi N_A \int_0^\infty \left[e^{-\beta u(r)} - 1\right] r^2 dr$$

| T (K) | $B_2$ (cm³/mol) Exp. | UA (TraPPE) | 5-sites |
|-------|---------------------|-------------|---------|
| 200 | –246 | –238 | –242 |
| 300 | –105 | –102 | –104 |
| 500 | –34 | –33 | –34 |

#### Propiedades de Transporte (Gas Diluido, 300 K)
| Propiedad | Exp. | MD (UA) | MD (5-sites) |
|-----------|------|---------|--------------|
| Viscosidad $\eta$ (μPa·s) | 11.2 | 11.0 | 11.1 |
| Difusividad $D$ (cm²/s) | 0.22 | 0.21 | 0.22 |
| Conductividad térmica $\lambda$ (mW/m·K) | 34.5 | 33.8 | 34.2 |

### 2. Fase Líquida

El metano líquido (T < 111.7 K a 1 atm) es un fluido molecular simple, **proto-tipo de líquido de Lennard-Jones**:

| Propiedad | Exp. (100 K) | MD TraPPE-UA | MD 5-sites |
|-----------|-------------|--------------|------------|
| Densidad (g/cm³) | 0.422 | 0.425 | 0.423 |
| $\Delta H_{vap}$ (kJ/mol) | 8.9 | 9.1 | 9.0 |
| Tensión superficial (mN/m) | 13.7 | 13.2 | 13.5 |
| Compresibilidad $\kappa_T$ (GPa⁻¹) | 2.1 | 2.0 | 2.1 |
| Viscosidad (μPa·s) | 118 | 110 | 115 |

#### Estructura Local: Función de Distribución Radial

$$g(r) = \frac{V}{N^2}\left\langle \sum_{i \neq j} \delta(r - r_{ij}) \right\rangle$$

- **Primer pico**: ~3.8 Å (contacto molecular)
- **Segundo pico**: ~7.0 Å (segundo shell)
- **Coordinación**: ~12–13 (empaquetamiento cercano, no tetraédrico)

A diferencia del agua, el metano líquido no presenta **estructura orientacional** significativa: las moléculas rotan casi libremente, y $g(r)$ depende solo de la distancia centro-centro.

### 3. Fase Sólida: Polimorfismo

El metano cristaliza en tres fases principales:

| Fase | Estructura | Rango de T | Grupo Espacial | Características |
|------|-----------|------------|----------------|-----------------|
| **I** | Cúbica (plástica) | > 20.4 K | Fm3m | Rotor libre, orientacionalmente desordenado |
| **II** | Hexagonal | 20.4–8.8 K | P6₃/mmc | Orden parcial, 8 moléculas por celda |
| **III** | Monoclínica | < 8.8 K | P2₁/c | Orden orientacional completo |

#### Fase Plástica (I)
- **Centros de masa**: red cúbica centrada en caras (fcc)
- **Orientación molecular**: desordenada, distribución isotrópica
- **Entropía de transición I→II**: ~6.0 J/mol·K
- **Dinámica**: reorientación molecular con tiempo de correlación ~1 ps

Las simulaciones MD con potenciales de orientación (Gay-Berne, Kihara) reproducen la fase plástica. El modelo UA puro **no captura** la transición de fase orientacional, requiriendo términos de anisotropía.

#### Fase III (Ordenada)
- Las moléculas adoptan orientaciones definidas
- Interacciones cuadrupolo-cuadrupolo ordenadas
- Barrera de reorientación: ~0.5–1.0 kcal/mol

---

## Espectroscopía Vibracional

### Modos Normales (3N – 6 = 9 modos vibracionales)

| Modo | Símbolo | Frecuencia (cm⁻¹) | Actividad | Descripción |
|------|---------|-------------------|-----------|-------------|
| $\nu_1$ | a₁ | 2917 | Raman | Simétrico C–H stretch |
| $\nu_2$ | e | 1534 | Raman | Degenerado H–C–H bend |
| $\nu_3$ | t₂ | 3019 | IR | Degenerado C–H stretch |
| $\nu_4$ | t₂ | 1306 | IR | Degenerado H–C–H bend |

### En Fase Condensada
| Efecto | Líquido (100 K) | Sólido I (20 K) |
|--------|----------------|-----------------|
| $\nu_1$ desplazamiento | –5 cm⁻¹ | –8 cm⁻¹ |
| $\nu_3$ ensanchamiento | ~15 cm⁻¹ | ~5 cm⁻¹ |
| $\nu_4$ actividad Raman | Aparece (ruptura de simetría) | Intensa |

La **ruptura de degeneración** en fase condensada (especialmente en fase III) indica distorsión del campo cristalino que reduce la simetría efectiva de Td a subgrupos.

---

## Mecánica Estadística del Metano

### Función de Partición
Para N moléculas en fase condensada:

$$Q = \frac{q_{trans}^{3N} \cdot q_{rot}^{3N} \cdot q_{vib}^{9N} \cdot q_{inter}^{N(N-1)/2}}{N! \Lambda^{3N}}$$

Donde:
- $q_{rot} = \frac{\sqrt{\pi}}{\sigma}\left(\frac{T^3}{\Theta_A \Theta_B \Theta_C}\right)^{1/2} = \frac{\sqrt{\pi}}{\sigma}\left(\frac{8\pi^2 I k_B T}{h^2}\right)^{3/2}$ (rotor esférico, $\sigma = 12$)
- $\Theta_{rot} = 7.54$ K (temperatura rotacional)

### Transición de Fase Líquido-Vapor
| Propiedad | Valor Crítico |
|-----------|--------------|
| $T_c$ | 190.56 K |
| $P_c$ | 4.599 MPa (45.99 bar) |
| $\rho_c$ | 0.1627 g/cm³ |
| $\omega$ (factor acéntrico) | 0.011 |

El metano tiene el **factor acéntrico más bajo** de todos los hidrocarburos, indicando comportamiento casi ideal. Las simulaciones GEMC (Gibbs Ensemble Monte Carlo) reproducen la curva de coexistencia con error < 2%.

---

## Hidratos de Metano: Estructura y Modelamiento

Los hidratos de gas son inclusiones cristalinas donde el metano se aloja en cavidades de una red de agua:

| Estructura | Cavidades | Radio (Å) | Ocupación típica | Condiciones de estabilidad |
|------------|-----------|-----------|-----------------|----------------------------|
| **sI** | 2 pequeñas (5¹²) + 6 grandes (5¹²6²) | 3.9 / 4.3 | 80–90% | P > 25 bar, T < 285 K |
| **sII** | 16 pequeñas (5¹²) + 8 grandes (5¹²6⁴) | 3.9 / 4.7 | 60–70% | P > 15 bar, T < 275 K |
| **sH** | 3 pequeñas + 2 medianas + 1 grande | — | Variable | Con ayudantes (C₃H₈, i-C₅H₁₂) |

### Modelamiento de Hidratos

**Potencial de interacción metano-agua**:
$$U_{CH_4-H_2O} = 4\epsilon_{CH_4-O}\left[\left(\frac{\sigma_{CH_4-O}}{r}\right)^{12} - \left(\frac{\sigma_{CH_4-O}}{r}\right)^{6}\right]$$

Con reglas de combinación:
$$\sigma_{CH_4-O} = \frac{\sigma_{CH_4} + \sigma_{O}}{2}, \quad \epsilon_{CH_4-O} = \sqrt{\epsilon_{CH_4} \cdot \epsilon_{O}}$$

**Simulaciones GCMC** predicen:
- Presión de disociación del hidrato sI: ~25 bar a 275 K (exp: ~26 bar)
- Calor de formación: ~54 kJ/mol H₂O
- Difusividad del metano en cavidades: ~10⁻⁹ cm²/s (extremadamente lenta)

---

## Adsorción y Almacenamiento

### 1. Carbón Activado y MOFs
- **Isoterma de Langmuir**: $\theta = \frac{KP}{1+KP}$
- **Isoterma de BET**: multicapa
- **Simulaciones GCMC**: predicen capacidad de adsorción a 298 K, 35 bar

| Material | Capacidad (mmol/g, 298 K, 35 bar) | Método |
|----------|----------------------------------|--------|
| Carbón activado | 8–12 | GCMC (UA) |
| MOF-5 | 15–18 | GCMC (AA) |
| ZIF-8 | 10–12 | GCMC (UA) |

### 2. Almacenamiento en Gas Natural (ANG)
- Densidad adsorbida en poros de 1–2 nm: ~200 g/L (comparable a gas comprimido a 250 bar)
- Simulaciones de DFT + MD predicen selectividad CH₄/CO₂ y CH₄/N₂

---

## Comparación con la Serie de Alcanos

| Propiedad | Metano | Etano | Butano |
|-----------|--------|-------|--------|
| **Simetría** | Td (esférica) | D₃d (cilíndrica) | C₂h/C₂v (flexible) |
| **Momento dipolar** | 0 | 0 | 0 |
| **Cuadrupolo** | Moderado | Alto | Muy alto |
| **Polarizabilidad** | Baja | Media | Alta |
| **Fase plástica** | Sí (I) | Sí (I) | No |
| **Barrera torsional** | N/A | ~2.9 kcal/mol | ~4.5 kcal/mol (gauche-anti) |
| **Modelo ideal** | UA (1 sitio) | UA (2 sitios) | AA o UA con torsión |
| **Comportamiento de fase** | Casi ideal ($\omega$=0.011) | Casi ideal | No ideal |
| **Temperatura crítica** | 190.6 K | 305.3 K | 425.1 K |

---

## Métodos Avanzados de Cálculo

### Dinámica *Ab Initio* para Metano
| Proceso | Método | Resultado |
|---------|--------|-----------|
| Disociación C–H | CCSD(T), DFT (PBE0) | Barrera ~105 kcal/mol |
| Metano + OH → CH₃ + H₂O | DFT (B3LYP) | Barrera ~6 kcal/mol |
| Metano en zeolita (adsorción) | DFT-D | Energía de adsorción ~15–25 kJ/mol |
| Metano + H₂O (hidrato) | MP2/aug-cc-pVTZ | Energía de estabilización ~12 kJ/mol |

### Machine Learning para Metano
- **SNAP, ACE, MACE**: potenciales equivariantes para simulaciones de metano en condiciones extremas (alta P, alta T)
- **Aplicación**: formación de hielos de metano en planetas gigantes (Neptuno, Urano)

---

## Resumen del Paradigma de Modelamiento

| Escala | Método | Propiedades |
|--------|--------|-------------|
| Electrónica (Å) | CCSD(T), DFT | Geometría, espectro, disociación |
| Molecular (nm) | MD/FF (TraPPE, OPLS) | EoS, transporte, fases |
| Mesoscópica (μm) | GCMC, DPD | Adsorción, hidratos, mezclas |
| Continuo (cm) | CFD + EoS (Peng-Robinson) | Ingeniería de procesos |

El metano es el **sistema de referencia** para:
1. **Parametrización de campos de fuerza UA**: si un modelo falla en metano, fallará en alcanos mayores
2. **Teoría de fluidos simples**: el metano líquido es el fluido de Lennard-Jones por excelencia
3. **Hidratos de gas**: sistema modelo para entender el confinamiento molecular
4. **Planetología**: hielos de metano en cuerpos del sistema solar exterior

> **Insight final**: la simplicidad estructural del metano —un átomo de carbono tetraédricamente rodeado de hidrógenos— esconde una riqueza fenomenológica en fase condensada. Su modelamiento exitoso requiere capturar no solo las interacciones de van der Waals, sino también las sutilezas electrostáticas (cuadrupolo) y, en sistemas confinados, los efectos de entropía configuracional.
