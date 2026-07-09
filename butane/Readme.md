 # Butano: Perspectiva desde el Modelamiento Molecular

## Estructura y Conformación

El butano (C₄H₁₀) es el alcano más simple que exhibe **conformacionalismo** significativo. Su cadena carbonada de cuatro átomos puede adoptar diversas conformaciones rotacionales alrededor del enlace C2–C3, cada una con distintas energías:

- **Anti** (o *trans*): conformación más estable. Los grupos metilo se encuentran a 180°. Mínimo de energía debido a la mínima repulsión estérica y torsional.
- **Gauche** (+60° y −60°): dos conformaciones degeneradas (mismo nivel energético). Los grupos metilo están a 60°, generando una repulsión estérica moderada (~0.9 kcal/mol por encima del anti).
- **Eclipsado** (sinclinal y anticlinal): máximos de energía. Los hidrógenos o grupos metilo se eclipsan, generando fuerte repulsión torsional.
- **Totalmente eclipsado** (0°): máximo absoluto de energía. Los dos grupos metilo se eclipsan directamente (~4.5–5.0 kcal/mol por encima del anti).

Esta distribución de energía conformacional se modela mediante potenciales dihedrales de la forma:

$$V(\phi) = \frac{V_n}{2}[1 + \cos(n\phi - \delta)]$$

Donde los parámetros $V_n$ (barreras rotacionales), $n$ (periodicidad) y $\delta$ (fase) se ajustan para reproducir la curva de energía potencial (PES) obtenida por cálculos *ab initio*.

---

## Propiedades Mecánicas desde el Modelamiento

### 1. **Dinámica Molecular (MD)**
En simulaciones de butano líquido o gas, las propiedades mecánicas emergen del promedio estadístico sobre el ensamble:

- **Presión y ecuación de estado**: se calcula mediante el virial de Clausius o fluctuaciones del volumen (NPT).
- **Viscosidad**: obtenida por el método de Green-Kubo (integrales de autocorrelación del tensor de presión) o ensayos de no-equilibrio (NEMD).
- **Coeficiente de difusión**: via la relación de Einstein o Green-Kubo.

### 2. **Mecánica Estadística**
La función de partición del butano debe incluir explícitamente las contribuciones conformacionales:

$$Q = q_{\text{trans}} \cdot q_{\text{rot}} \cdot q_{\text{vib}} \cdot q_{\text{conf}}$$

Donde $q_{\text{conf}}$ suma sobre los estados anti y gauche con sus respectivos factores de Boltzmann. Esto afecta directamente:
- Capacidad calorífica ($C_p$, $C_v$)
- Entropía configuracional
- Cambio de entalpía en transiciones de fase

### 3. **Propiedades Elásticas en Fase Condensada**
En butano líquido a alta presión (relevante para aplicaciones como propelente o combustible):
- **Módulo de compresibilidad** $K_T = -V(\partial P/\partial V)_T$
- **Coeficiente de expansión térmica** $\alpha_P$

Estos se extraen de simulaciones NPT o de la pendiente de la isoterma en ecuaciones de estado (como PC-SAFT o Peng-Robinson).

---

## Parámetros de Campo de Fuerzas (Force Field)

Para modelar el butano clásicamente, los campos de fuerza como **OPLS-AA**, **CHARMM** o **AMBER** definen:

| Parámetro | Valor típico (OPLS-AA) | Significado físico |
|-----------|------------------------|-------------------|
| $r_0$ (C–C) | 1.53 Å | Longitud de enlace |
| $r_0$ (C–H) | 1.09 Å | Longitud de enlace |
| $k_b$ (C–C) | ~268 kcal/mol·Å² | Constante de rigidez |
| $k_\theta$ (H–C–H) | ~37 kcal/mol·rad² | Rigidez angular |
| $V_3$ (torsión C–C) | ~2.0 kcal/mol | Barrera rotacional |
| $\sigma_{CC}$ | 3.50 Å | Diámetro de Lennard-Jones |
| $\epsilon_{CC}$ | 0.066 kcal/mol | Pozo de potencial LJ |

Los parámetros de **no-enlace** (Lennard-Jones 12-6) gobiernan las propiedades mecánicas colectivas:
$$U_{LJ}(r) = 4\epsilon\left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right]$$

---

## Aplicaciones del Modelo

| Propiedad | Método de simulación | Resultado esperado |
|-----------|---------------------|-------------------|
| Densidad del líquido (298 K, 1 atm) | MD NPT | ~0.573 g/cm³ |
| Presión de vapor | MD/GCMC o EoS | ~2.1 atm a 298 K |
| Viscosidad | NEMD o Green-Kubo | ~0.17 cP (gas, 298 K) |
| Tensión superficial | MD con método de la gota | ~12 mN/m |

---

## Perspectiva Cuántica vs. Clásica

- **Cálculos DFT/HF** (B3LYP, MP2): determinan con precisión la curva de PES conformacional, las barreras rotacionales y las frecuencias vibracionales. Son la base para **parametrizar** los campos de fuerza clásicos.
- **Dinámica *ab initio*** (CPMD, BOMD): útil para estudiar la ruptura de enlaces o reacciones (combustión, pirolisis), donde los efectos electrónicos dominan.

---

## Resumen

Desde el modelamiento molecular, el butano es un sistema **paradigmático** para entender:
1. **Conformacionalismo**: la interconversión anti/gauche como fundamento de la termodinámica configuracional.
2. **Transferibilidad de parámetros**: los grupos metilo y metileno del butano se "transfieren" a alcanos mayores.
3. **Escala de tiempo**: las rotaciones alrededor de C–C ocurren en ~ps, mientras que la difusión molecular en líquido ocurre en ~ns.
