# Este script extrae la energía total y los parámetros de celda de un archivo de salida
# de un cálculo de DFT. La energía total se utiliza para calcular la energía de
# formación de un material. Los parámetros de celda se utilizan para calcular
# las propiedades cristalinas de un material. La energía total se obtiene a partir
# de la energía de formación del material. Los parámetros de celda se obtienen
# a partir de la estructura cristalina del material. La energía de formación se
# puede utilizar para calcular la estabilidad de un material. La energía de
# formación se puede utilizar para calcular la energía de un material a partir
# de la energía de formación de los elementos que lo componen. 

import re
import numpy as np

def parse_nwchem_output(file_path):
    energy = None
    lattice_params = []
    
    with open(file_path, 'r') as f:
        for line in f:
            # Extraer energía total
            if "Total DFT energy" in line:
                energy = float(re.search(r'-\d+\.\d+', line).group())
            
            # Extraer parámetros de celda
            if "Lattice parameters" in line:
                params = re.findall(r'\d+\.\d+', line)
                lattice_params = list(map(float, params))
    
    return energy, lattice_params

# Ejemplo de uso
energy, lattice = parse_nwchem_output("optimizacion.out")
print(f"Energia total: {energy} eV")
print(f"Parámetros de celda: {lattice}")