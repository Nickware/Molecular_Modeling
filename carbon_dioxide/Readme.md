 # Dióxido de Carbono (CO₂): Perspectiva desde el Modelamiento Molecular

## Estructura y Geometría

El CO₂ es una molécula **lineal** con simetría D∞h, lo que simplifica enormemente su descripción mecánica pero introduce complejidades en sus propiedades colectivas:

- **Distancia C=O**: ~1.16 Å
- **Ángulo O–C–O**: 180° (estrictamente lineal)
- **Momento dipolar**: cero (simetría centro-inversión)

Esta linealidad impone restricciones cinemáticas en dinámica molecular: solo dos grados de libertad rotacionales (no tres) y un modo de vibración asimétrico activo en IR.

---

## Descripción Cuántica y Orbitales

Desde DFT o Hartree-Fock:

- **Orbitales**: el carbono presenta hibridación sp, con dos orbitales σ formados por solapamiento sp–p de oxígeno, y dos enlaces π deslocalizados perpendiculares al eje molecular.
- **Carga formal**: C(+4), O(–2), pero el análisis de población (Mulliken, NBO, Bader) revela una polarización parcial: C ≈ +1.2, O ≈ –0.6 en campos de fuerza clásicos.

La linealidad se mantiene incluso con correlación electrónica (CCSD(T)), aunque en fase condensada las interacciones intermoleculares pueden inducir **ligeras distorsiones** (< 1°) en cristales bajo alta presión.

---

## Campo de Fuerzas para CO₂

### Modelo de Sitios Múltiples (TraPPE, EPM2)

Dado que CO₂ no tiene momento dipolar permanente, los campos de fuerza **átomo-átomo** puros fallan en reproducir sus propiedades de fase. Los modelos exitosos incluyen **sitios de carga fuera de los núcleos**:

| Parámetro | TraPPE | EPM2 | Significado |
|-----------|--------|------|-------------|
| $q_C$ | +0.6512 e | +0.6512 e | Carga parcial en C |
| $q_O$ | –0.3256 e | –0.3256 e | Carga parcial en O |
| $q_{M}$ (sitio central) | — | –0.6512 e | Carga en sitio ficticio |
| $\sigma_{OO}$ | 3.05 Å | 3.033 Å | Diámetro LJ oxígeno |
| $\epsilon_{OO}$ | 79.0 K | 80.507 K | Pozo LJ oxígeno |
| $l_{C-M}$ | — | 0.699 Å | Distancia C–sitio carga |

El **modelo de 3 sitios de carga** (EPM2) coloca una carga negativa en el centro de la molécula para reproducir el **cuadrupolo** ($\Theta = -4.3 \times 10^{-26}$ esu·cm²), esencial para las interacciones intermoleculares.

### Potencial Intermolecular

$$U_{ij} = 4\epsilon_{ij}\left[\left(\frac{\sigma_{ij}}{r_{ij}}\right)^{12} - \left(\frac{\sigma_{ij}}{r_{ij}}\right)^{6}\right] + \frac{q_i q_j}{4\pi\epsilon_0 r_{ij}}$$

Con reglas de combinación de Lorentz-Berthelot para los parámetros cruzados.

---

## Propiedades Mecánicas desde el Modelamiento

### 1. **Fase Gaseosa**
- **Ecuación de estado virial**: los coeficientes $B_2(T)$, $B_3(T)$ se calculan por integración del potencial de par:
$$B_2(T) = -2\pi N_A \int_0^\infty [e^{-\beta u(r)} - 1] r^2 dr$$
- **Viscosidad y difusión**: Chapman-Enskog con potencial de Stockmayer (LJ + cuadrupolo).

### 2. **Fase Líquida y Supercrítica**
El CO₂ supercrítico es un **fluido modelo** en la industria:

| Propiedad | Método MD | Resultado (304.25 K, 73.8 bar) |
|-----------|-----------|-------------------------------|
| Densidad | NPT | ~0.468 g/cm³ (punto crítico) |
| Viscosidad | Green-Kubo | ~0.03 cP |
| Difusividad | Einstein | ~10⁻⁴ cm²/s |
| Tensión superficial | Gota/Placa | ~0 (supercrítico) |

La **compresibilidad isotérmica** $\kappa_T$ diverge en el punto crítico, un desafío numérico que requiere tamaños de sistema grandes y técnicas de remuestreo.

### 3. **Fase Sólida (Hielo Seco)**
El CO₂ cristaliza en estructura **Pa3** (cúbica simple, grupo espacial 205):

- **Sublimación**: a 1 atm, sublima a 194.7 K. La entalpía de sublimación (~25.2 kJ/mol) se reproduce bien con modelos de 3 sitios.
- **Módulo de bulk** $K_0$: ~8.8 GPa (experimental). Simulaciones NPT con potenciales flexibles lo predicen dentro del 10%.
- **Transiciones de fase de alta presión**: a > 12 GPa, CO₂ molecular transforma a fases poliméricas (fase III, IV, V), requiriendo **dinámica *ab initio*** (DFT-PBE, BLYP) para describir la polimerización C=O → C–O–C.

---

## Propiedades Espectroscópicas y Vibracionales

| Modo | Frecuencia (exp.) | Actividad | Descripción |
|------|-------------------|-----------|-------------|
| $\nu_1$ (simétrico) | 1333 cm⁻¹ | Raman | Estiramiento C=O simétrico |
| $\nu_2$ (deformación) | 667 cm⁻¹ | IR | Flexión degenerada (πᵤ) |
| $\nu_3$ (asimétrico) | 2349 cm⁻¹ | IR | Estiramiento C=O asimétrico |

En simulaciones MD con espectroscopía de correlación:
$$I(\omega) \propto \int_{-\infty}^{\infty} \langle \mu(0) \cdot \mu(t) \rangle e^{i\omega t} dt$$

El ensanchamiento de líneas en fase condensada (líquida/sólida) revela **dephasing** por colisiones y reorientación molecular.

---

## CO₂ como Fluido de Referencia en Termodinámica

### Ecuaciones de Estado (EoS)
- **Peng-Robinson**: parámetro $a(T)$ con función alfa de Soave; $\omega_{CO_2} = 0.225$.
- **Span-Wagner**: EoS multiparamétrica con 42 coeficientes, referencia IAPWS para CO₂ industrial.

### Modelos de Agregación
El CO₂ forma **clatratos** (hidratos) con agua a alta presión. Las simulaciones de Monte Carlo en ensamble gran-canónico (GCMC) predicen las condiciones de estabilidad de las estructuras sI y sII.

---

## Mecánica Cuántica Avanzada

### Cálculos *Ab Initio*
- **MP2/aug-cc-pVTZ**: reproduce la geometría y frecuencias vibracionales con error < 1%.
- **CCSD(T)-F12**: "gold standard" para la energía de disociación $D_0$ = 532 kJ/mol (promedio de los dos enlaces C=O).
- **DFT**: funcionales híbridos (B3LYP, PBE0) son adecuados para geometría; funcionales con dispersión (B97-D, ωB97X-D) necesarios para interacciones intermoleculares.

### Dinámica *Ab Initio* (CPMD, BOMD)
- **Disociación**: a T > 3000 K, CO₂ se disocia → CO + O. Las barreras se calculan por NEB (Nudged Elastic Band).
- **Reacción con H₂O**: formación de ácido carbónico H₂CO₃, un intermediario fugaz en la captura de CO₂.

---

## Resumen del Paradigma de Modelamiento

| Escala | Método | Propiedades |
|--------|--------|-------------|
| Electrónica | HF, DFT, CCSD(T) | Estructura, espectro, energías de disociación |
| Molecular | FF (EPM2, TraPPE) | EoS, transporte, fases |
| Mesoscópica | DPD, LBM | Flujo en medios porosos (captura geológica) |
| Continuo | CFD con EoS | Inyección supercrítica, almacenamiento |

El CO₂ es un **sistema de prueba fundamental**: su simplicidad estructural contrasta con la riqueza de sus propiedades colectivas (punto crítico, supercríticidad, polimerización a alta presión), haciéndolo ideal para validar metodologías de modelamiento molecular.