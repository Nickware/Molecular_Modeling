#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Análisis de Convergencia de Energía DFT de NWChem
=============================================

Descripción:
-----------
Este programa extrae y visualiza datos de convergencia de energía DFT de archivos de salida de NWChem.
Grafica la energía SCF versus el número de iteraciones, mostrando el comportamiento de convergencia.

Características:
-----------------
- Analiza archivos de salida de NWChem para valores "Total DFT energy"
- Convierte energías de Hartree a electronvolts (eV)
- Genera gráficos de convergencia
- Anota el valor final de energía convergida

Uso:
-----

1. Coloque archivo de salida NWChem (por defecto: 'optimizacion.nwo') en la misma carpeta
2. Ejecute script: python dft_convergence_plotter.py
3. Salida: 'dft_energy_convergencia.png' gráfico

Autor: N.Torres
Fecha: 2025-04-27
Versión: 1.0
Licencia: MIT
"""


import matplotlib.pyplot as plt
import numpy as np
import re

# Función para extraer energías DFT de un archivo NWChem
def extract_energies(filename):
    """Extrae Energías DFT de archivo NWChem"""
    energies = []
    with open(filename, 'r') as f:
        for line in f:
            if "Total DFT energy" in line:
                # Handle various whitespace and formatting possibilities
                match = re.search(r"-\d+\.\d+", line)
                if match:
                    try:
                        energy = float(match.group())
                        energies.append(energy)
                    except ValueError:
                        continue
    return energies

# Parametros
output_file = "optimizacion.nwo"  # Cambie al archivo de salida NWChem de su elección
hartree_to_ev = 27.2114       # Factor de conversión

# Extrae y procesa energías
energies_ha = extract_energies(output_file)
if not energies_ha:
    raise ValueError("No se encontraron energías DFT en el archivo de salida!")

# Convierte a eV y calcula números de iteraciones
energies_ev = [e * hartree_to_ev for e in energies_ha]
iterations = list(range(1, len(energies_ev) + 1))

# Crea figura
plt.figure(figsize=(10, 6))

# Gráfico principal de energía
plt.plot(iterations, energies_ev, 'bo-', markersize=6, linewidth=1.5, label="Energía DFT")
plt.xlabel("Número de iteración SCF", fontsize=12)
plt.ylabel("Energía (eV)", fontsize=12)
plt.grid(True, linestyle='--', alpha=0.5)

# Anotaaciones
converged_iteration = iterations[-1]
converged_energy = energies_ev[-1]
plt.annotate(f'Converged: {converged_energy:.6f} eV\nat iteration {converged_iteration}',
             xy=(converged_iteration, converged_energy),
             xytext=(10, 10), textcoords='offset points',
             bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
             arrowprops=dict(arrowstyle='->'))

# Title and legend
plt.title("Convergencia de Energía DFT\nOptimización SCF NWChem", fontsize=14, pad=20)
plt.legend(loc='upper right')

plt.tight_layout()
plt.savefig("dft_energy_convergence.png", dpi=300, bbox_inches="tight")
plt.show()