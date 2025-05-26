# =========================================================
# Script: graficar_optimizacion_geométrica.py
# Autor: N.Torres
# Descripción: 
#   Este script lee un archivo de salida de NWChem, extrae las energías
#   de cada iteración de una optimización geométrica, convierte las
#   energías de Hartree a electronvoltios (eV) y grafica la convergencia.
# Fecha de creación: 2025-04-27
# =========================================================

import re
import matplotlib.pyplot as plt

# Funciones básicas para el script

def leer_archivo(ruta):
    with open(ruta, 'r') as f:
        return f.read()

def extraer_energias(contenido):
    # Buscar energías en cada paso de optimización
    energias = []
    for linea in contenido.splitlines():
        if "Total DFT energy" in linea:
            match = re.search(r'Total DFT energy\s+=\s+([-\d\.E]+)', linea)
            if match:
                energia_hartree = float(match.group(1))
                energias.append(energia_hartree)
    return energias

def hartree_a_ev(energia_hartree):
    # 1 Hartree = 27.2114 eV
    return energia_hartree * 27.2114

def graficar_optimizacion(energias_hartree):
    iteraciones = list(range(1, len(energias_hartree) + 1))
    energias_ev = [hartree_a_ev(e) for e in energias_hartree]

    plt.figure(figsize=(8, 5))
    plt.plot(iteraciones, energias_ev, marker='o')
    plt.xlabel('Iteración')
    plt.ylabel('Energía Total (eV)')
    plt.title('Convergencia de la Optimización Geométrica')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# Programa principal
if __name__ == "__main__":
    # Cambiar a la ruta que se requiera
    ruta_archivo = 'optimizacion.nwo'  
    contenido = leer_archivo(ruta_archivo)
    energias = extraer_energias(contenido)

    if energias:
        graficar_optimizacion(energias)
    else:
        print("No se encontraron energías de optimización en el archivo.")
