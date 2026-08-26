# Structure optimization of butane molecule with NWChem calculator
# (Simplified version without detailed input parameters)
# Note: Ensure NWChem is properly installed and configured in your environment.
# You may need to adjust the calculator settings based on your NWChem installation.

import numpy as np
import matplotlib.pyplot as plt

from ase import Atoms
from ase.io import read, write
from ase.optimize import BFGS
from ase.calculators.nwchem import NWChem
from ase.vibrations import Vibrations

# ---------------------------------------------------------------------------
# 1. Optimización de geometría
# ---------------------------------------------------------------------------
butane = read('butane.xyz')
butane.calc = NWChem(xc='PBE')

opt = BFGS(butane, trajectory="optimization.traj")
opt.run(fmax=0.01)

write('butane_opt.xyz', butane)

energy = butane.get_potential_energy()
print(f"Energía potencial en el mínimo: {energy:.6f} eV")

# ---------------------------------------------------------------------------
# 2. Fuerzas atómicas en la geometría optimizada
# ---------------------------------------------------------------------------
forces = butane.get_forces()
max_force = np.max(np.linalg.norm(forces, axis=1))
print(f"Fuerza máxima residual sobre un átomo: {max_force:.6f} eV/Å")

# ---------------------------------------------------------------------------
# 3. Análisis vibracional
# ---------------------------------------------------------------------------
vib = Vibrations(butane)
vib.run()
vib.summary()

# Frecuencias en cm^-1 (ASE devuelve complejos: la parte imaginaria
# codifica los modos con curvatura negativa, es decir, frecuencias
# "imaginarias" en la convención espectroscópica habitual)
frequencies = vib.get_frequencies()  # array complejo, unidades cm^-1

real_part = frequencies.real
imag_part = frequencies.imag

# Una frecuencia es "imaginaria" si su parte imaginaria es
# significativamente distinta de cero (dominante sobre la parte real)
is_imaginary = imag_part > 1e-6
n_imaginary = np.sum(is_imaginary)

# Para graficar, representamos los modos imaginarios como valores
# negativos en el eje (convención estándar en espectros vibracionales)
plot_values = np.where(is_imaginary, -imag_part, real_part)

# ---------------------------------------------------------------------------
# 4. Gráfica diagnóstica: ¿mínimo verdadero o punto de silla?
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9, 5))

colors = ['crimson' if imag else 'steelblue' for imag in is_imaginary]
indices = np.arange(1, len(plot_values) + 1)

ax.bar(indices, plot_values, color=colors, edgecolor='black', linewidth=0.5)
ax.axhline(0, color='black', linewidth=0.8)

ax.set_xlabel("Modo vibracional")
ax.set_ylabel(r"Frecuencia (cm$^{-1}$)")
ax.set_title("Espectro vibracional — diagnóstico de mínimo de energía")
ax.set_xticks(indices)

# Leyenda manual
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='steelblue', edgecolor='black', label='Frecuencia real (curvatura positiva)'),
    Patch(facecolor='crimson', edgecolor='black', label='Frecuencia imaginaria (curvatura negativa)')
]
ax.legend(handles=legend_elements, loc='best')

# Anotación de conclusión directamente en la gráfica
if n_imaginary == 0:
    veredicto = "✓ Mínimo de energía verdadero (0 frecuencias imaginarias)"
    box_color = "#d4edda"
elif n_imaginary == 1:
    veredicto = "⚠ Punto de silla de primer orden (1 frecuencia imaginaria)"
    box_color = "#fff3cd"
else:
    veredicto = f"✗ Punto de silla de orden {n_imaginary} ({n_imaginary} frecuencias imaginarias)"
    box_color = "#f8d7da"

ax.text(
    0.5, -0.22, veredicto,
    transform=ax.transAxes,
    ha='center', va='top',
    fontsize=11, fontweight='bold',
    bbox=dict(boxstyle="round,pad=0.4", facecolor=box_color, edgecolor='gray')
)

plt.tight_layout()
plt.savefig("diagnostico_frecuencias.png", dpi=200, bbox_inches='tight')
plt.show()

# ---------------------------------------------------------------------------
# 5. Resumen en consola
# ---------------------------------------------------------------------------
print("\n--- Diagnóstico de la geometría optimizada ---")
print(f"Número de frecuencias imaginarias: {n_imaginary}")
if n_imaginary == 0:
    print("La estructura corresponde a un mínimo de energía verdadero.")
else:
    print("La estructura NO es un mínimo verdadero: corresponde a un punto")
    print("de silla. Revise la geometría inicial o el camino de optimización.")
