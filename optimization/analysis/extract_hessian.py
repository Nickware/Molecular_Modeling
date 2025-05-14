# Este script extrae la matriz Hessiana de un archivo de salida de una optimización
# de energía de un cálculo de DFT. La matriz Hessiana se utiliza para calcular
# las constantes elásticas de un material. La matriz Hessiana se obtiene a partir
# de la segunda derivada de la energía con respecto a las coordenadas atómicas.
# La matriz Hessiana se puede utilizar para calcular las constantes elásticas
# de un material a partir de la relación entre la energía y las coordenadas atómicas.

import numpy as np

def read_hessian(file_path):
    hessian = []
    in_hessian = False
    
    with open(file_path, 'r') as f:
        for line in f:
            if "Hessian of energy" in line:
                in_hessian = True
                continue
            if in_hessian and line.strip() == "":
                break
            if in_hessian:
                hessian.extend(list(map(float, line.strip().split())))
    
    size = int(len(hessian)**0.5)
    return np.array(hessian).reshape(size, size)

# Ejemplo de cálculo de constantes elásticas
hessian = read_hessian("optimizacion.out")
volume = 125.6  # Volumen en Å³
Cij = hessian * (1e21) / (volume * 1e-30)  # Conversión a Pa