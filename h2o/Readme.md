# Agua (H₂O): Perspectiva desde el Modelamiento Molecular

## Estructura y Geometría de Equilibrio

El agua es la molécula **más estudiada en la historia de la química computacional**, no por su simplicidad, sino por la extraordinaria complejidad que emerge de ella.

### Geometría de Gas Aislado
| Parámetro | Valor Experimental | CCSD(T)/CBS |
|-----------|-------------------|-------------|
| Distancia O–H | 0.9572 Å | 0.958 Å |
| Ángulo H–O–H | 104.52° | 104.5° |
| Momento dipolar $\mu$ | 1.855 D | 1.85 D |

La geometría angular (no lineal) y la alta electronegatividad del oxígeno generan:
- **Momento dipolar permanente**: ~1.85 D (dominado por los pares solitarios del oxígeno)
- **Polarizabilidad**: $\alpha = 1.47 \times 10^{-30}$ m³ (anisotrópica)
- **Cuadrupolo**: $\Theta_{zz} = -2.63$ DA

---

## El Desafío del Enlace de Hidrógeno

El agua no se modela como molécula aislada. Su comportamiento está gobernado por el **puente de hidrógeno** (HB), una interacción intermolecular con características parcialmente covalentes:

### Geometría del Puente de Hidrógeno
| Parámetro | Valor Típico |
|-----------|-------------|
| Distancia O···O | 2.75–2.85 Å (líquido), 2.76 Å (hielo Ih) |
| Distancia O–H···O | ~1.8 Å (líquido) |
| Ángulo O–H···O | ~170° (casi lineal) |
| Energía de HB | ~5–6 kcal/mol (cooperativo: hasta 8–10 kcal/mol en cadenas) |

### Naturaleza del Puente de Hidrógeno
El análisis de población (NBO, ELF, QTAIM) revela:
- **Transferencia de carga**: ~0.03–0.05 e del aceptor al dador
- **Interacción orbital**: donación n(O) → σ*(O–H)
- **Componente electrostática**: ~50% de la energía total
- **Dispersión y correlación**: ~25% (crítica, requiere métodos con dispersión)
- **Intercambio-repulsión**: ~25%

> **Insight clave**: el HB no es puramente electrostático. La correlación electrónica (dispersión) y la transferencia de carga son esenciales. Los modelos puramente puntiformes (cargas fijas) fallan en capturar la **cooperatividad** y la **fluctuación**.

---

## Modelos de Campo de Fuerza para Agua

### 1. Modelos de Cargas Fijas (Fixed-Charge)

| Modelo | Año | Carga O (e) | $\sigma_{OO}$ (Å) | $\epsilon_{OO}$ (kcal/mol) | Características |
|--------|-----|-------------|-------------------|---------------------------|---------------|
| **SPC** | 1981 | –0.82 | 3.166 | 0.155 | 3 sitios, rígida |
| **SPC/E** | 1987 | –0.8476 | 3.166 | 0.155 | Corrección polarización media |
| **TIP3P** | 1983 | –0.834 | 3.151 | 0.152 | Optimizado para biofísica |
| **TIP4P** | 1983 | — | 3.154 | 0.155 | 4 sitios (carga M en bisectriz) |
| **TIP4P/2005** | 2005 | — | 3.159 | 0.185 | Reproduce propiedades de fase |
| **TIP5P** | 2000 | — | 3.120 | 0.160 | 5 sitios (tetraédrico) |

**Ecuación de potencial (TIP4P)**:
$$U = \sum_{i<j} \frac{q_i q_j}{r_{ij}} + 4\epsilon\left[\left(\frac{\sigma}{r_{OO}}\right)^{12} - \left(\frac{\sigma}{r_{OO}}\right)^{6}\right]$$

### 2. Modelos Polarizables

| Modelo | Tipo | Características |
|--------|------|---------------|
| **POL3** | Fluctuating charge | Cargas que responden al entorno |
| **SPC-FQ** | Fluctuating charge | Carga O: –0.82 → –0.95 en fase condensada |
| **AMOEBA** | Induced dipole | Polarizabilidad atómica explícita |
| **iAMOEBA** | Induced dipole + multipolos | Hasta cuadrupolos inducidos |
| **MB-pol** | Many-body expansion | Hasta términos de 3 cuerpos |

**Energía de polarización**:
$$U_{pol} = -\frac{1}{2}\sum_i \boldsymbol{\mu}_i^{ind} \cdot \mathbf{E}_i^{0}$$

Donde $\boldsymbol{\mu}_i^{ind} = \alpha_i \mathbf{E}_i^{tot}$ y el campo total incluye contribuciones de todas las moléculas.

### 3. Modelos *Ab Initio* (Machine Learning)

| Modelo | Base | Precisión |
|--------|------|-----------|
| **GAP** (Gaussian Approximation Potential) | DFT (PBE0-TS) | Energía ~1 meV/H₂O |
| **NNP** (Neural Network Potential) | CCSD(T) | Energía ~0.5 meV/H₂O |
| **MACE** | Equivariante E(3) | Transferible, rápido |
| **q-TIP4P/F** | Path-integral + reparametrización | Efectos cuánticos nucleares |

---

## Propiedades Mecánicas desde el Modelamiento

### 1. Fase Líquida: El Agua "Anómala"

El agua líquida exhibe propiedades que desafían los fluidos normales:

| Propiedad | Experimental | SPC/E | TIP4P/2005 | MB-pol |
|-----------|-----------|-------|-----------|--------|
| Densidad (298 K, 1 atm) | 0.997 g/cm³ | 0.998 | 0.999 | 0.997 |
| $\Delta H_{vap}$ | 10.5 kcal/mol | 11.2 | 10.6 | 10.5 |
| $\epsilon$ (constante dieléctrica) | 78.4 | 71 | 60 | 78 |
| $\kappa_T$ (GPa⁻¹) | 0.46 | 0.49 | 0.44 | 0.46 |
| $D$ (difusividad, 10⁻⁵ cm²/s) | 2.3 | 2.5 | 1.8 | 2.2 |
| $\eta$ (viscosidad, mPa·s) | 0.89 | 0.5 | 0.85 | 0.88 |

#### Estructura Local: Funciones de Distribución Radial

$$g_{OO}(r), \quad g_{OH}(r), \quad g_{HH}(r)$$

- **$g_{OO}(r)$**: primer pico a 2.75 Å (tetraédrico), segundo pico a 4.5 Å (segundo shell tetraédrico)
- **$g_{OH}(r)$**: pico intramolecular a 0.96 Å, pico intermolecular a 1.8 Å (HB)
- **Número de coordinación**: ~4.5 (no tetraédrico estricto, ~80% tetraédrico a 298 K)

#### Red de Puentes de Hidrógeno
- **Fracción de HB**: ~70–80% en agua líquida ambiente
- **Vida media de un HB**: ~1–2 ps
- **Reorganización cooperativa**: la ruptura de un HB induce reordenamiento de 3–4 moléculas vecinas

#### Anomalías Mecánicas
| Anomalía | Origen Molecular | Modelo que la Reproduce |
|----------|---------------|----------------------|
| **Máximo de densidad a 4°C** | Compromiso entre empaquetamiento tetraédrico (baja densidad) y empaquetamiento cercano (alta densidad) | TIP4P/2005, MB-pol |
| **Expansión al congelar** | Hielo Ih: estructura tetraédrica abierta, $\rho = 0.917$ g/cm³ | Todos los modelos tetraédricos |
| **Aumento de viscosidad con T** (a baja T) | Aumento de conectividad de la red de HB | Modelos polarizables |
| **Disminución de $\kappa_T$ con T** (0–46°C) | Transición de estructura "baja densidad" (LDL) a "alta densidad" (HDL) | TIP4P/2005, modelos de 2 estados |

### 2. Fase Sólida: Polimorfismo del Hielo

El agua tiene **19 fases cristalinas** conocidas (más amorfos), cada una con distinta topología de red de HB:

| Fase | Estructura | Densidad (g/cm³) | Condiciones | Grupo Espacial |
|------|-----------|------------------|-------------|----------------|
| **Ih** | Hexagonal | 0.917 | 1 atm, T < 273 K | P6₃/mmc |
| **Ic** | Cúbica | 0.931 | Metaestable | Fd3m |
| **II** | Romboédrica | 1.17 | > 0.2 GPa | R3 |
| **III** | Tetragonal | 1.16 | 0.2–0.4 GPa | P4₁2₁2 |
| **V** | Monoclínica | 1.23 | 0.4–0.6 GPa | A2/a |
| **VI** | Tetragonal | 1.31 | 0.6–2.1 GPa | P4₂/nmc |
| **VII** | Cúbica (body-centered) | 1.50 | > 2.1 GPa | Pn3m |
| **VIII** | Tetragonal | 1.51 | > 2.1 GPa, T < 273 K | I4₁/amd |
| **X** | Cúbica (ionica) | ~2.5 | > 70 GPa | Pn3m |

#### Hielo Ih: Estructura Tetraédrica
- Cada O coordinado tetraédricamente por 4 O vecinos
- Los H están desordenados (regla de Bernal-Fowler): cada O tiene 2 H cercanos (covalentes) y 2 H lejanos (HB)
- Entropía residual: $R\ln(3/2) \approx 3.4$ J/mol·K (Pauling)

#### Transiciones de Fase
- **Ih → II**: reconstrucción de red, colapso parcial de la estructura abierta
- **VI → VII/VIII**: transición de red molecular a red con puente de H simétrico (O–H–O lineal, H en el centro)
- **VII → X**: ionización completa, O²⁻ con H⁺ en sitios intersticiales (superiónico)

### 3. Fase Supercrítica

| Propiedad | Valor Crítico | Comportamiento |
|-----------|--------------|----------------|
| $T_c$ | 647.1 K | |
| $P_c$ | 22.06 MPa | |
| $\rho_c$ | 0.322 g/cm³ | |
| $\kappa_T$ | Diverge | Fluctuaciones gigantes de densidad |

En el régimen supercrítico, el agua pierde la red de HB extendida pero mantiene **agregados locales** de 3–5 moléculas. Las simulaciones MD muestran:
- Fracción de HB: ~20–30% a 700 K, 30 MPa
- Difusividad: ~10⁻³ cm²/s (similar a líquido a baja presión)

---

## Dinámica Molecular y Efectos Cuánticos Nucleares

### El Problema del Protón

El hidrógeno es ligero ($m_H \approx 1$ u). A 300 K, su longitud de onda de de Broglie térmica:

$$\lambda_{dB} = \frac{h}{\sqrt{2\pi m_H k_B T}} \approx 1.5 \text{ Å}$$

Esto es **comparable** a la distancia O–H (0.96 Å). Los efectos cuánticos nucleares (NQE) son **no despreciables**:

| Efecto | Magnitud | Método de Simulación |
|--------|----------|----------------------|
| **Túnel del protón** | Barreras < 5 kcal/mol: significativo | Path-Integral MD (PIMD) |
| **Cero-punto de energía** | ~5.5 kcal/mol por modo O–H | PIMD, RPMD |
| **Efecto isotópico** | D₂O: $\rho_{max}$ a 11.2°C vs 4°C para H₂O | PIMD con masas isotópicas |
| **Difusividad** | D₂O: ~20% más lento que H₂O | RPMD, CMD |

### Path-Integral Molecular Dynamics (PIMD)

El protón se representa como un **polímero de cuentas** (beads) conectadas por resortes:

$$H_{PI} = \sum_{i=1}^{P} \left[\frac{p_i^2}{2m^*} + \frac{1}{2}m\omega_P^2(x_i - x_{i+1})^2 + \frac{1}{P}V(x_i)\right]$$

Donde $\omega_P = P/\beta\hbar$ y $P$ es el número de beads (típicamente 32–64 a 300 K).

**Resultados clave de PIMD**:
- El protón en el puente de H tiene una distribución de probabilidad **más ancha** que la clásica
- El HB se **debilita** por efectos de cero-punto (el protón "fluctúa" más)
- La densidad del agua líquida a 300 K es ~2% mayor con NQE que clásicamente

---

## Espectroscopía Vibracional

### Modos Normales de la Molécula Aislada

| Modo | Frecuencia (cm⁻¹) | Actividad | Descripción |
|------|-------------------|-----------|-------------|
| $\nu_1$ (a₁) | 3657 | Raman | Simétrico O–H stretch |
| $\nu_2$ (a₁) | 1595 | IR/Raman | Bend H–O–H |
| $\nu_3$ (b₂) | 3756 | IR | Asimétrico O–H stretch |

### En Fase Condensada

| Región | Frecuencia (cm⁻¹) | Asignación |
|--------|-------------------|------------|
| **Stretch O–H** | 3200–3600 | Bandas de HB: red desordenada |
| **Bend H–O–H** | 1640 | Intramolecular, poco afectado |
| **Libración** | 400–1000 | Rotación inhibida (rocking, wagging, twisting) |
| **Stretch O···O** | 150–200 | Modo de red intermolecular |
| **Bend O–O–O** | < 100 | Modo acústico |

#### Ensanchamiento y Desplazamiento
- **Desplazamiento al rojo**: $\nu_{OH}$ disminuye de 3657 cm⁻¹ (gas) a ~3400 cm⁻¹ (líquido) debido al debilitamiento del enlace O–H por donación de carga en el HB
- **Ensanchamiento**: FWHM ~300–400 cm⁻¹ en líquido (vs ~10 cm⁻¹ en gas), por inhomogeneidad del entorno y dephasing ultrarrápido (~50 fs)

---

## Propiedades Mecánicas Detalladas

### 1. Tensión Superficial

$$\gamma = \frac{1}{2A}\left\langle \sum_{i<j} \left(r_{ij} - \frac{3z_{ij}^2}{r_{ij}}\right) \frac{dU}{dr_{ij}} \right\rangle$$

| T (K) | $\gamma$ (mN/m) | Modelo MD |
|-------|-----------------|-----------|
| 273 | 75.6 | TIP4P/2005: 74.2 |
| 298 | 71.7 | TIP4P/2005: 70.5 |
| 373 | 58.9 | TIP4P/2005: 57.8 |

El perfil de densidad en la interfaz líquido-vapor muestra una **capa de enriquecimiento** molecular (~3 Å de espesor) con orientación preferencial: dipolos paralelos a la superficie.

### 2. Viscosidad

La viscosidad del agua es anómalamente alta para una molécula tan pequeña:

$$\eta = \frac{V}{k_B T} \int_0^\infty \langle P_{xy}(0)P_{xy}(t) \rangle dt$$

Donde $P_{xy}$ es el componente off-diagonal del tensor de presión.

| T (K) | $\eta$ (mPa·s) | Mecanismo Molecular |
|-------|---------------|---------------------|
| 273 | 1.79 | Red de HB rígida, cooperativa |
| 298 | 0.89 | Ruptura/reformación de HB |
| 373 | 0.28 | Red fragmentada, comportamiento más "normal" |

Los modelos polarizables (AMOEBA, iAMOEBA) reproducen $\eta$ mucho mejor que los de carga fija, indicando que la **respuesta electrónica** es clave para el transporte de momentum.

### 3. Difusividad

$$D = \frac{1}{6}\lim_{t\to\infty} \frac{d}{dt}\langle |\mathbf{r}(t) - \mathbf{r}(0)|^2 \rangle$$

| T (K) | $D$ (10⁻⁵ cm²/s) | Mecanismo |
|-------|-----------------|-----------|
| 273 | 1.1 | Saltos entre pozos de HB |
| 298 | 2.3 | Saltos + difusión continua |
| 350 | 6.0 | Difusión predominante |

El mecanismo de difusión no es un "salto" discreto sino una **reorganización cooperativa**: una molécula rompe 2 HB, rota ~60°, y forma nuevos HB con vecinos.

---

## Agua en Sistemas Confinados e Interfaciales

### 1. Hidratación de Biomoléculas

| Entorno | Comportamiento del Agua | Modelo Requerido |
|---------|------------------------|----------------|
| **Primera capa de solvatación** (proteínas) | Tiempo de residencia ~10–100 ps, orientación preferencial | Polarizable, explícito |
| **Groove del ADN** | Puentes de H con bases, estabilización estructural | Polarizable, explícito |
| **Interior de micelas** | Agua "moleculas de bulto" con propiedades alteradas | UA + agua explícita |
| **Canales de membrana** | Filtración por tamaño, selectividad iónica | Polarizable, PIMD |

### 2. Agua en Nanoporos

- **Carbono**: hidrofobicidad extrema, contact angle > 90°
- **Sílice**: hidrofilicidad, capa de agua adsorbida ~3–5 Å
- **Zeolitas**: confinamiento unidimensional, difusividad anisotrópica

Las simulaciones de **non-equilibrium MD (NEMD)** predicen flujos que desafían la ley de Hagen-Poiseuille debido a la **deslización en la pared** (slip length).

---

## Métodos de Cálculo *Ab Initio* para Agua

### DFT: El Problema de la Dispersión

| Funcional | Descripción Dispersión | Precisión en Agua |
|-----------|----------------------|-------------------|
| LDA | No incluye | HB sobreestimado, estructura incorrecta |
| PBE | No incluye | Densidad líquida ~15% alta |
| BLYP | No incluye | Similar a PBE |
| PBE-D3 | Grimme D3 | Mejora significativa |
| B97-D | Incluida intrínsecamente | Buena para clusters |
| vdW-DF | Funcional no local | Densidad líquida correcta |
| SCAN | Meta-GGA | Bueno sin corrección explícita |
| **r²SCAN** | Meta-GGA mejorado | Excelente para fases |

### Métodos de Referencia

| Método | Precisión | Costo | Aplicación |
|--------|-----------|-------|------------|
| MP2 | ~1 kcal/mol en HB | $O(N^5)$ | Clusters hasta 20 H₂O |
| CCSD(T) | ~0.1 kcal/mol | $O(N^7)$ | Clusters hasta 10 H₂O |
| DLPNO-CCSD(T) | ~0.2 kcal/mol | $O(N)$ aprox. | Clusters grandes |
| Quantum Monte Carlo | Exacta (en principio) | Muy alto | Benchmarks |

### Dinámica *Ab Initio* (CPMD, BOMD)

- **CPMD**: Car-Parrinello, electrones ficticios. Time-step ~0.1 fs.
- **BOMD**: Born-Oppenheimer, electrones convergidos. Time-step ~0.5 fs.
- **Tiempo de simulación típico**: 10–100 ps (vs 100 ns para FF clásicos)

**Aplicaciones clave**:
- Disociación autoionización: 2 H₂O → H₃O⁺ + OH⁻ (barrera ~20 kcal/mol en bulk)
- Transferencia de protón en cadenas de HB (Grotthuss mechanism)
- Reacciones en superficies (electrocatálisis, fotocatálisis)

---

## Resumen: La Jerarquía del Modelamiento del Agua

| Nivel | Método | Propiedades | Limitaciones |
|-------|--------|-------------|--------------|
| **Clásico rígido** | SPC/E, TIP3P | Estructura, termodinámica aproximada | No NQE, no polarización, viscosidad incorrecta |
| **Clásico flexible** | SPC/Fw, q-TIP4P/F | Espectro, NQE aproximado | Cargas fijas |
| **Polarizable** | AMOEBA, iAMOEBA | Dieléctrica, viscosidad, solvatación | Costo ~10×, parametrización compleja |
| **Machine Learning** | MB-pol, MACE | Energía CCSD(T), NQE | Costo ~100–1000×, transferibilidad limitada |
| **Ab initio MD** | DFT-D, SCAN | Química reactiva, proton transfer | Tiempo de simulación corto, dispersión |
| **PIMD + ab initio** | RPMD, PIMD@DFT | NQE exactos, isótopos | Costo extremo (~10⁴×) |

---

## El Agua como Sistema Paradigmático

El agua encapsula los **desafíos fundamentales** del modelamiento molecular:

1. **Interacciones no aditivas**: el potencial de 2 cuerpos no basta (cooperatividad de HB)
2. **Polarización dinámica**: las cargas no son fijas
3. **Efectos cuánticos nucleares**: el protón no es una partícula clásica
4. **Múltiples escalas**: desde el femtosegundo de la vibración O–H hasta el segundo de la difusión en poros

Cada avance en la modelización del agua —desde los primeros modelos de Bernal y Fowler (1933) hasta los potenciales de machine learning de hoy— ha sido un **test de fuego** para la química computacional.

> *"El agua es la sustancia más estudiada y menos entendida."* — atribuido a diversos autores, y aún vigente.