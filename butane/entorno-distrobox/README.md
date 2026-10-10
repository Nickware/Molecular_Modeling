# Ejecución del ejemplo de butano en Distrobox

Esta guía configura un contenedor Ubuntu 24.04 con Distrobox para ejecutar `butane.py`. El repositorio permanece en el directorio personal del host, compartido con el contenedor por Distrobox; el entorno virtual de Python se guarda aparte del repositorio.

## Requisitos en el host

- Linux con Distrobox y un proveedor de contenedores compatible, como Podman o Docker.
- Git, si todavía no tienes una copia local del repositorio.
- Conexión a internet para descargar la imagen y los paquetes.
- Espacio y memoria suficientes para los cálculos DFT. El tiempo depende del equipo y de la instalación de NWChem.

Instala Distrobox y Podman con el gestor de paquetes de tu distribución y comprueba que estén disponibles:

```bash
distrobox --version
podman --version
```

Docker puede usarse en lugar de Podman si Distrobox está configurado para ese proveedor.

## Crear y entrar al contenedor

```bash
distrobox create --name butane-dft --image ubuntu:24.04
distrobox enter butane-dft
```

Los pasos siguientes se ejecutan dentro del contenedor.

## Instalar NWChem y las dependencias de Python

```bash
sudo apt update
sudo apt install -y nwchem python3 python3-venv python3-pip
```

Verifica que NWChem esté en el `PATH`:

```bash
command -v nwchem
apt-cache policy nwchem
```

Si `apt` no encuentra el paquete `nwchem`, comprueba que el componente `universe` de Ubuntu esté habilitado, actualiza el índice de paquetes y vuelve a intentar la instalación.

Crea un entorno virtual e instala ASE, NumPy y Matplotlib:

```bash
python3 -m venv ~/.venvs/butane
source ~/.venvs/butane/bin/activate
python -m pip install --upgrade pip
python -m pip install ase numpy matplotlib
```

## Ejecutar el cálculo

Desde el contenedor, cambia al directorio `butane` del repositorio. Ajusta la ruta si clonaste el proyecto en otra ubicación:

```bash
cd ~/gitHub/MolecularModeling/butane
source ~/.venvs/butane/bin/activate
MPLBACKEND=Agg python3 butane.py
```

`MPLBACKEND=Agg` selecciona un backend gráfico sin interfaz, adecuado para ejecutar dentro del contenedor incluso si no hay una sesión gráfica disponible. La gráfica se guarda igualmente como `diagnostico_frecuencias.png`.

El script lee `butane.xyz`, ejecuta la optimización geométrica con BFGS y NWChem, calcula las frecuencias vibracionales e imprime energía, fuerza residual y diagnóstico. Al finalizar, revisa `butane_opt.xyz`, `optimization.traj` y `diagnostico_frecuencias.png` en esa misma carpeta. ASE puede dejar además archivos de caché vibracional.

## Problemas frecuentes

- **`nwchem` no se encuentra:** confirma que está instalado dentro del contenedor y que `command -v nwchem` muestra su ruta.
- **No se puede instalar `nwchem`:** habilita `universe` en los repositorios de Ubuntu y ejecuta de nuevo `sudo apt update`.
- **No se encuentra `butane.xyz`:** ejecuta el programa desde `butane/`; el script usa una ruta relativa al directorio actual.
- **Error al importar un módulo de Python:** activa el entorno virtual con `source ~/.venvs/butane/bin/activate` y repite la instalación de paquetes.
- **La ejecución tarda o no converge:** revisa la salida de NWChem y ASE. El tiempo y la convergencia dependen de los recursos disponibles y de los parámetros electrónicos, que en el script se dejan parcialmente a los valores predeterminados.

Para conocer qué calcula cada etapa y qué no está cubierto por este ejemplo, consulta el [README del proyecto de butano](../Readme.md).