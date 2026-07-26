 # Heptano (C₇H₁₆): Perspectiva desde el Modelamiento Molecular

## Estructura e Isomería

El heptano marca la transición hacia los **alcanos de cadena media**, donde la complejidad conformacional y las propiedades colectivas comienzan a exhibir comportamientos que anticipan los de polímeros y lubricantes.

### Isómeros Estructurales

El C₇H₁₆ posee **9 isómeros estructurales**:

| Isómero | Nomenclatura IUPAC | Simetría | Punto de Ebullición (K) | Punto de Fusión (K) |
|---------|-------------------|----------|------------------------|---------------------|
| **n-heptano** | Heptano | C₂ (conformación all-anti) | 371.6 | 182.6 |
| **2-metilhexano** | Isoheptano | C₁ | 363.2 | 154.9 |
| **3-metilhexano** | — | C₁ | 365.0 | — |
| **2,2-dimetilpentano** | Neoheptano (uno de ellos) | C₂v | 352.6 | 149.4 |
| **2,3-dimetilpentano** | — | C₁ | 362.9 | — |
| **2,4-dimetilpentano** | — | C₂v | 353.7 | 154.0 |
| **3,3-dimetilpentano** | — | C₂v | 359.7 | 138.7 |
| **2,2,3-trimetilbutano** | — | C₃v | 354.0 | 248.6 |
| **3-etilpentano** | — | C₂ | 366.6 | — |

> **Insight clave**: el n-heptano es el **fluido de referencia** para el índice de cetano (definido como 100), lo que lo convierte en el estándar molecular para modelar la **ignición por compresión** en motores de combustión interna.

---

## n-Heptano: Conformacionalismo en 4D

### Grados de Libertad Torsionales

El n-heptano posee **cuatro ángulos dihedrales internos** rotacionales (C2–C3, C3–C4, C4–C5, C5–C6), generando un **paisaje energético de 4 dimensiones** de extraordinaria riqueza.

### Enumeración de Rotámeros

Si cada dihedral adopta 3 estados (t, g⁺, g⁻), existen **3⁴ = 81 conformaciones** teóricas. Sin embargo, las restricciones estéricas reducen drásticamente las poblables:

| Tipo de Conformación | Notación | Energía Relativa (kcal/mol) | Población (~300 K) |
|---------------------|----------|----------------------------|-------------------|
| **all-trans** (tttt) | C₂v (aprox) | 0.00 | ~35% |
| **Un gauche** (tgtt, ttgt, etc.) | — | ~0.5 | ~40% (suma) |
| **Dos gauche adyacentes opuestos** (tg⁺g⁻t) | — | ~1.0 | ~12% |
| **Dos gauche adyacentes sinclinales** (tg⁺g⁺t) | — | ~3.5–4.0 | < 0.5% |
| **Tres o más gauche** | — | ~1.5–5.0 | < 5% |

### Reglas de Rotámero y Cooperatividad

Para cadenas de ≥ 5 carbonos, emergen **reglas empíricas** que los campos de fuerza deben reproducir:

1. **Regla del gauche-pentano**: dos giros gauche adyacentes en la misma dirección (g⁺g⁺ o g⁻g⁻) generan una penalización de ~3.5 kcal/mol por interacción 1,5.

2. **Regla del gauche-hexano**: en cadenas de 6+ carbonos, las interacciones 1,6 comienzan a ser relevantes. Para heptano, la secuencia g⁺tg⁺ en posiciones 2-3-4-5 genera un contacto C1···C6 desfavorable.

3. **Efecto cooperativo**: la población de ttgt no es exactamente el producto de las poblaciones individuales de t y g debido a correlaciones conformacionales a lo largo de la cadena.

### Longitud de Persistencia y Flexibilidad

En el límite de cadena larga, la longitud de persistencia caracteriza la rigidez:

$$l_p = \frac{\langle \mathbf{r}_i \cdot \mathbf{r}_{i+1} \rangle}{1 - \cos\theta}$$

Para el n-heptano líquido a 298 K:
- **Longitud de persistencia**: $l_p \approx$ 7.5–8.0 Å
- **Radio de giro**: $R_g \approx$ 3.4–3.6 Å
- **Longitud de contorno total** (all-trans): ~9.2 Å

Esto sitúa al heptano en la frontera entre **molécula rígida** y **cadena flexible**: aún conserva memoria orientacional, pero ya exhibe plegamiento significativo.

---

## Geometría de Equilibrio y Parámetros Estructurales

### Conformación all-trans (Mínimo Global)

| Parámetro | Valor (Exp./ED) | B3LYP/6-31G** | MP2/cc-pVTZ |
|-----------|----------------|---------------|-------------|
| r(C–C) | 1.534 Å | 1.531 Å | 1.533 Å |
| r(C–H) | 1.096 Å | 1.093 Å | 1.094 Å |
| ∠(C–C–C) | 113.2° | 113.0° | 113.1° |
| ∠(H–C–H) | 107.8° | 107.5° | 107.7° |
| Dihedral C–C–C–C | 180.0° | 179.8° | 179.9° |

### Distorsiones Conformacionales

En la conformación **tg⁺tt**:
- El ángulo C–C–C adyacente al gauche se **abre** a ~114.5° para aliviar la repulsión 1,4
- Los enlaces C–H orientados hacia el interior del pliegue se acortan ~0.002 Å por tensión estérica

---

## Campo de Fuerzas para Heptano

### Modelos de Átomos Unidos (UA)

El heptano es el **sistema de prueba estándar** para validar campos de fuerza UA en hidrocarburos de cadena media:

| Modelo | $\sigma_{CH_2}$ (Å) | $\epsilon_{CH_2}$ (K) | $\sigma_{CH_3}$ (Å) | $\epsilon_{CH_3}$ (K) | Torsión $V_3$ (kcal/mol) |
|--------|---------------------|----------------------|---------------------|----------------------|--------------------------|
| **TraPPE-UA** | 3.93 | 46.0 | 3.75 | 98.0 | 2.0 (transferido) |
| **NERD** | 3.93 | 47.0 | 3.75 | 100.0 | 2.0 |
| **SKS** | 3.905 | 47.5 | 3.73 | 102.0 | 2.1 |
| **MARTINI (CG)** | 4.7 (bead) | 120 (bead) | — | — | Potencial de ángulo |

### Modelos All-Atom (AA)

| Modelo | Carga C (e) | Carga H (e) | Torsión Cruzada | 1,5-VDW |
|--------|-------------|-------------|-----------------|---------|
| **OPLS-AA** | –0.12 | +0.06 | No | Escala 1.0 |
| **CHARMM** | –0.18 | +0.06 | No | Escala 1.0 |
| **AMBER (GAFF)** | –0.12 | +0.06 | No | Escala 1.0 |
| **MM4** | –0.15 | +0.075 | Sí (explícita) | Escala 1.0 |

### Desafíos de Parametrización Específicos del Heptano

1. **Términos de torsión cruzada**: para 4 dihedrales, se requieren términos de acoplamiento $\phi_i \times \phi_j$ para reproducir la estadística conformacional correcta.

2. **Interacciones 1,5 y 1,6**: en la conformación all-trans, la distancia C1···C5 es ~4.0 Å (sin contacto). En conformaciones plegadas, C1···C6 puede acercarse a ~3.5 Å, requiriendo tratamiento explícito de no-enlace.

3. **Transferibilidad**: los parámetros optimizados para butano/pentano deben predecir heptano sin reajuste. TraPPE-UA logra esto con error < 2% en densidad.

---

## Propiedades Mecánicas desde el Modelamiento

### 1. Fase Líquida

El heptano líquido es el **fluido modelo** para gasolina, disolventes orgánicos y lubricantes base.

| Propiedad | Exp. (298 K, 1 atm) | MD TraPPE-UA | MD OPLS-AA |
|-----------|---------------------|--------------|------------|
| Densidad (g/cm³) | 0.684 | 0.683 | 0.689 |
| $\Delta H_{vap}$ (kJ/mol) | 36.6 | 36.2 | 37.0 |
| Tensión superficial (mN/m) | 20.1 | 19.5 | 20.8 |
| Viscosidad (mPa·s) | 0.397 | 0.35 | 0.42 |
| Difusividad (10⁻⁵ cm²/s) | 3.2 | 3.5 | 2.9 |
| Compresibilidad $\kappa_T$ (GPa⁻¹) | 1.45 | 1.40 | 1.48 |
| Calor específico $C_p$ (J/mol·K) | 224.8 | 218 | 230 |

#### Estructura Local y Orden Parcial

$$g_{CC}(r), \quad g_{CH_2-CH_2}(r)$$

- **Primer pico $g_{CC}$**: ~5.8 Å (mayor que pentano: ~5.2 Å)
- **Coordinación de primeros vecinos**: ~9–10 (menor que metano por forma alargada)
- **Orden orientacional**: el parámetro de orden $Q_2 \approx 0.20$ indica **dominios nemáticos locales** (cybotactic clusters) de ~15–25 Å

> A diferencia de alcanos cortos, el heptano líquido desarrolla una **estructura de capas** (layering) en la interfaz y en confinamiento, con cadenas alineadas paralelamente.

#### Mecanismo de Difusión

$$D = \frac{1}{6} \lim_{t\to\infty} \frac{d}{dt} \langle |\mathbf{r}(t) - \mathbf{r}(0)|^2 \rangle$$

El heptano exhibe **reptación** (reptation) como mecanismo dominante:
- Tiempo de residencia en un "tubo" de vecinos: ~10–20 ps
- Saltos conformacionales (trans-gauche): ~1–2 ps
- Reorientación end-to-end: ~100 ps

### 2. Fase Gaseosa

#### Segundo Coeficiente Virial

| T (K) | $B_2$ (cm³/mol) Exp. | TraPPE-UA | OPLS-AA |
|-------|---------------------|-----------|---------|
| 350 | –1230 | –1190 | –1250 |
| 400 | –780 | –760 | –790 |
| 500 | –380 | –370 | –385 |

La forma alargada del heptano reduce drásticamente $B_2$, indicando fuertes interacciones intermoleculares efectivas incluso a moderadas presiones.

#### Viscosidad del Gas

A 400 K, 1 atm: $\eta \approx 7.5$ μPa·s. Las simulaciones de MD en gas diluido con potenciales UA reproducen este valor con error < 5%.

### 3. Fase Sólida y Cristalización

El n-heptano cristaliza en fase **ortorrómbica**:

| Parámetro | Valor |
|-----------|-------|
| Grupo espacial | Pna2₁ o Pbca (dependiendo de polimorfo) |
| Moléculas por celda | 4 |
| Conformación en sólido | all-trans plano-zigzag |
| Densidad (sólido, 150 K) | ~0.82 g/cm³ |
| $\Delta H_{fus}$ | 14.0 kJ/mol |

#### Desafío del Subenfriamiento en MD

La cristalización del heptano es **lenta** (cinéticamente controlada). En simulaciones de MD:
- El líquido se **subenfría** fácilmente hasta 100 K por debajo del punto de fusión
- La nucleación homogénea requiere tiempos de microsegundos (inaccesibles en MD atómico)
- Técnicas de **metadinámica** o **umbrella sampling** se usan para estudiar la barrera de nucleación

---

## Espectroscopía Vibracional

### Modos Normales del n-Heptano (Conformación all-trans, C₂h)

Con 23 átomos: 3N – 6 = **63 modos vibracionales**:

| Región (cm⁻¹) | N° de Modos | Actividad | Descripción |
|---------------|-------------|-----------|-------------|
| 2800–3000 | 18 | Raman + IR | Stretching C–H simétrico/asimétrico |
| 1400–1500 | 12 | Raman + IR | Deformación C–H (scissoring, rocking) |
| 1300–1400 | 8 | Raman + IR | Wagging C–H |
| 1000–1200 | 8 | Raman + IR | Twisting + stretching C–C |
| 800–1000 | 6 | Raman + IR | Rocking + modos de cadena |
| 400–600 | 4 | Raman | Modos de cadena (chain modes) |
| 100–300 | 4 | Raman | Torsión C–C (libración interna) |
| < 100 | 3 | Raman | Modos de red (traslación, libración) |

### Modos de Cadena (Chain Modes)

Característicos de alcanos de cadena media/larga:

| Modo | Frecuencia (cm⁻¹) | Descripción |
|------|-------------------|-------------|
| **Longitudinal Acoustic** (LA) | ~200–400 | Compresión/expansión de la cadena completa |
| **Twisting Acoustic** (TA) | ~100–200 | Torsión colectiva de la cadena |
| **Bending Acoustic** (BA) | ~50–150 | Flexión de la cadena |

Estos modos colectivos son **sensibles a la longitud de cadena** y se desplazan a menor frecuencia conforme aumenta el número de carbonos.

---

## Comparación con la Serie de Alcanos Lineales

| Propiedad | Pentano | Hexano | Heptano | Octano |
|-----------|---------|--------|---------|--------|
| **N° de rotores** | 2 | 3 | 4 | 5 |
| **Dimensiones PES** | 2D | 3D | **4D** | 5D |
| **Población all-trans** | ~55% | ~48% | ~35% | ~28% |
| **Longitud de persistencia** | ~5.5 Å | ~6.5 Å | **~7.5 Å** | ~8.5 Å |
| **$R_g$ (líquido)** | ~2.6 Å | ~3.0 Å | **~3.5 Å** | ~3.9 Å |
| **Viscosidad (mPa·s, 298K)** | 0.22 | 0.30 | **0.40** | 0.54 |
| **Tensión superficial (mN/m)** | 15.5 | 17.9 | **20.1** | 21.8 |
| **$T_c$ (K)** | 469.7 | 507.6 | **540.2** | 568.7 |
| **$P_c$ (bar)** | 33.7 | 30.3 | **27.4** | 24.9 |
| **Factor acéntrico $\omega$** | 0.252 | 0.296 | **0.350** | 0.398 |

> **Tendencia**: a medida que aumenta la longitud de cadena, el alcano se vuelve más **flexible** (menor población all-trans), más **viscoso**, y más **no ideal** (mayor factor acéntrico).

---

## Aplicaciones Industriales y Modelamiento

### 1. Índice de Cetano y Combustión

El n-heptano es el **estándar de ignición** (cetano = 100). Desde el modelamiento molecular:

- **Dinámica de Reacción Directa (DRS)**: estudio de la disociación C–C y C–H en condiciones de motor
  - Barrera C–C: ~85 kcal/mol
  - Barrera C–H (abstracción por OH): ~6 kcal/mol
- **Mecanismos cinéticos**: simulaciones de ReaxFF para la oxidación en fase gas
- **Atomización**: simulaciones de MD con potenciales reactivos predicen la formación de radicales CH₃ y C₂H₅

### 2. Lubricación y Fluidos Base

El heptano modela los **grupos alquilo** en aceites lubricantes:
- **Película delgada**: simulaciones de MD no-equilibrio (NEMD) predicen viscosidad bajo cizallamiento
- **Confinamiento**: en películas de < 10 nm entre superficies, el heptano exhibe **solidificación por confinamiento** (aumento de viscosidad > 10×)

### 3. Adsorción y Separación

| Material | Aplicación | Método de Simulación |
|----------|-----------|----------------------|
| Zeolitas 5A | Separación n/i-parafinas | GCMC + MD |
| MOF-5 | Almacenamiento | GCMC |
| Carbón activado | Purificación | GCMC, DFT-D |
| Sílice mesoporosa | Cromatografía | MD, DFT |

Las simulaciones de **configuracional-bias Monte Carlo (CBMC)** son esenciales para muestrear las conformaciones del heptano en poros estrechos (~1 nm), donde las cadenas largas deben **plegarse** para entrar.

### 4. Coarse-Graining (MARTINI, SAFT-γ Mie)

Para simulaciones de mezclas heptano/agua o heptano/asfaltenos:
- **MARTINI**: 3 beads por heptano (1 terminal + 2 centrales)
- **SAFT-γ Mie**: modelo de cadenas tangentes de esferas Mie con parámetros derivados de termodinámica

---

## Métodos Avanzados de Cálculo

### Dinámica *Ab Initio* para Heptano

| Propiedad | Método | Resultado |
|-----------|--------|-----------|
| Barrera anti-gauche | MP2/cc-pVTZ | ~0.5 kcal/mol |
| Penalización g⁺g⁺ (1,5) | CCSD(T)/CBS | ~3.4 kcal/mol |
| Penalización g⁺tg⁺ (1,6) | DFT (B3LYP) | ~1.2 kcal/mol |
| Modos de cadena (LA) | DFT (ωB97X-D) | ~350 cm⁻¹ |

### Machine Learning para Heptano

- **ANI-2x**: predice energías de conformaciones con error ~0.3 kcal/mol, permitiendo exploración del PES 4D
- **MACE, NequIP**: potenciales equivariantes E(3) para simulaciones de heptano en condiciones extremas (alta P, alta T)
- **TorchMD**: simulaciones de 1000+ moléculas de heptano con precisión DFT a velocidad de FF clásico

---

## Resumen del Paradigma de Modelamiento

| Escala | Método | Propiedades | Desafíos Específicos del Heptano |
|--------|--------|-------------|----------------------------------|
| Electrónica (Å) | DFT, MP2, CCSD(T) | Geometría, barreras, espectro | PES 4D, cooperatividad, 1,6-interacciones |
| Molecular (nm) | MD/FF (TraPPE, OPLS) | EoS, transporte, fases | Cristalización lenta, viscosidad, orden local |
| Mesoscópica (μm) | DPD, MesoDyn, MARTINI | Mezclas, microemulsiones, asfaltenos | CG de cadenas semiflexibles |
| Continuo (cm) | CFD + EoS (PC-SAFT, Peng-Robinson) | Motores, separaciones, lubricación | Isomería, ramificación, cetano |

El heptano representa el **umbral de la complejidad conformacional** en alcanos:
- Con **4 rotores torsionales**, el espacio conformacional ya no es explorable por fuerza bruta sin métodos inteligentes (metadinámica, CBMC)
- Con **longitud de persistencia comparable a su tamaño**, exhibe **orden local nemático** en fase líquida
- Como **estándar de cetano**, conecta el modelamiento molecular con la **ingeniería de combustión**

> **Insight final**: el heptano es donde la química orgánica computacional deja de ser "molécula pequeña" y comienza a ser "cadena polimérica". Su modelamiento exitoso requiere no solo capturar la química de enlaces y torsiones, sino también la física estadística de cadenas flexibles, la termodinámica de fases condensadas, y la cinética de procesos de transporte.