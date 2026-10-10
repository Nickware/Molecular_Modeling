# Instalación del entorno ASE y sus calculadores

Esta guía configura un entorno Linux con Distrobox para los ejemplos del repositorio que usan [ASE](https://ase-lib.org/). Instala ASE y el backend de química cuántica NWChem, y comprueba que ambos pueden comunicarse con un cálculo pequeño de H2.

ASE es la interfaz Python; no sustituye al programa que calcula energías y fuerzas. La mayoría de los ejemplos ASE de este repositorio usan NWChem. El ejemplo `carbon_dioxide/quantum_expresso` usa Quantum ESPRESSO, que requiere instalar y configurar ese backend y sus pseudopotenciales por separado; esta guía no los instala.

## 1. Preparar el host

Distrobox necesita Linux y un motor de contenedores. Se recomienda Podman en modo rootless; Docker también es compatible si el usuario puede ejecutarlo sin `sudo`.

Ejecuta en una terminal del sistema anfitrión:

```bash
cat /etc/os-release
```

Si aparece un runtime como `Freedesktop SDK` en lugar de la distribución Linux (por ejemplo, en el terminal de VS Code instalado como Flatpak), instala y administra Distrobox desde una terminal del host, fuera del sandbox de la aplicación.

Instala Distrobox y Podman con el gestor de paquetes de tu distribución. Ejemplos:

```bash
# Ubuntu 22.10 o posterior
sudo apt update
sudo apt install -y distrobox podman

# Fedora
sudo dnf install -y distrobox podman

# Arch Linux
sudo pacman -Syu distrobox podman
```

Para otras distribuciones, consulta la [guía oficial de instalación y compatibilidad de Distrobox](https://distrobox.it/installation/). Comprueba que el motor funciona como usuario normal:

```bash
distrobox --version
podman --version
podman info
```

`podman info` debe terminar sin errores. Con Docker, comprueba `docker info` y configura los permisos siguiendo su [documentación oficial](https://docs.docker.com/engine/install/linux-postinstall/). Distrobox integra el contenedor con el host y comparte el directorio personal; no es un sandbox de seguridad.

## 2. Crear el contenedor

El repositorio debe estar dentro de tu directorio personal para que Distrobox lo comparta con el contenedor. Si aún no tienes una copia local, clónala en el host:

```bash
mkdir -p ~/gitHub
cd ~/gitHub
git clone https://github.com/Nickware/Molecular_Modeling.git MolecularModeling
```

Si ya tienes el repositorio en otro directorio bajo `$HOME`, no hace falta clonarlo de nuevo. Crea el contenedor en el host y entra:

```bash
distrobox create --name molecular-modeling --image ubuntu:24.04
distrobox enter --name molecular-modeling
```

La primera creación descarga la imagen y puede tardar. Continúa los siguientes pasos dentro de la sesión Ubuntu del contenedor.

## 3. Instalar ASE y NWChem

Comprueba que puedes instalar paquetes dentro del contenedor:

```bash
whoami
sudo -v
```

Escoge una sola de las siguientes opciones para instalar ASE, sus dependencias y NWChem.

### Opción A: Conda-forge

Miniforge proporciona Conda configurado para conda-forge. El instalador se selecciona según la arquitectura del contenedor:

```bash
sudo apt update
sudo apt install -y ca-certificates curl
curl -L "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$(uname -m).sh" -o /tmp/Miniforge3.sh
bash /tmp/Miniforge3.sh -b -p "$HOME/miniforge3"
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda create -y --override-channels -c conda-forge -n ase-nwchem python=3.12 ase nwchem numpy matplotlib
conda activate ase-nwchem
```

En una sesión futura, vuelve a cargar `conda.sh` y activa `ase-nwchem` antes de usar Python o NWChem.

### Opción B: repositorios de Ubuntu 24.04

En Ubuntu 24.04, NWChem y ASE están en el componente `universe`. Instala además las bibliotecas Python que requieren los ejemplos:

```bash
sudo apt update
sudo apt install -y nwchem python3-ase python3-numpy python3-matplotlib
```

Comprueba que APT encontró ambos paquetes:

```bash
apt-cache policy nwchem python3-ase
```

Si no aparece una versión candidata, habilita `universe` y repite la instalación:

```bash
sudo apt install -y software-properties-common
sudo add-apt-repository -y universe
sudo apt update
sudo apt install -y nwchem python3-ase python3-numpy python3-matplotlib
```

Con esta opción usa `python3` del sistema; los módulos instalados por APT podrían no estar disponibles desde un entorno virtual separado.

Si `sudo` no existe dentro del contenedor o el usuario no tiene privilegios, no ejecutes `sudo distrobox`. Consulta `distrobox enter --help` para el modo root admitido por el motor configurado o revisa la inicialización de la imagen.

## 4. Configurar ASE para NWChem

ASE necesita localizar el ejecutable `nwchem` para que su calculador lo lance. Confirma que está instalado y visible en el entorno activo:

```bash
command -v nwchem
```

Configura el comando de NWChem que usará ASE. `PREFIX` se sustituye por el nombre base de los archivos de entrada y salida:

```bash
export ASE_NWCHEM_COMMAND="nwchem PREFIX.nwi > PREFIX.nwo"
```

Los scripts basados en este backend usan `ase.calculators.nwchem.NWChem`; por ejemplo, `NWChem(xc="PBE")`. Esta variable de entorno debe definirse en cada sesión antes de ejecutar los ejemplos. Para confirmar que ASE importa la interfaz y puede crear un calculador:

```bash
python3 - <<'PY'
from ase import __version__
from ase.calculators.nwchem import NWChem

print("ASE:", __version__)
print("Calculador:", NWChem(xc="PBE", basis="3-21G"))
PY
```

## 5. Prueba de instalación ASE-NWChem

Ejecuta un cálculo de energía sencillo para H2. Este test no depende de ninguno de los proyectos moleculares; verifica Python, ASE, el comando configurado y el ejecutable NWChem juntos. Los archivos que NWChem genere quedan en un directorio temporal:

```bash
python3 - <<'PY'
import os
import tempfile

from ase import Atoms
from ase.calculators.nwchem import NWChem

with tempfile.TemporaryDirectory(prefix="ase-nwchem-test-") as workdir:
    os.chdir(workdir)
    molecule = Atoms("H2", positions=[(0, 0, 0), (0, 0, 0.74)])
    molecule.calc = NWChem(xc="PBE", basis="3-21G")
    energy = molecule.get_potential_energy()
    print(f"Prueba ASE-NWChem correcta: energía H2 = {energy:.6f} eV")
PY
```

El test debe terminar sin traceback e imprimir una energía numérica. Si falla, comprueba que activaste el entorno elegido, que `command -v nwchem` funciona y que `ASE_NWCHEM_COMMAND` está exportada en esa misma sesión. No continúes con un cálculo más grande hasta resolver el error.

## 6. Ejecutar un ejemplo del repositorio

Entra en la carpeta del proyecto que vas a ejecutar. Por ejemplo, para butano:

```bash
cd ~/gitHub/MolecularModeling/butane
```

Para cualquier proyecto que use NWChem, prepara el entorno elegido y ejecuta su script desde la carpeta que contiene sus archivos de entrada. Con Conda:

```bash
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate ase-nwchem
export ASE_NWCHEM_COMMAND="nwchem PREFIX.nwi > PREFIX.nwo"
MPLBACKEND=Agg python3 butane.py
```

Con APT, omite la activación de Conda y ejecuta el mismo comando de configuración y Python del sistema:

```bash
export ASE_NWCHEM_COMMAND="nwchem PREFIX.nwi > PREFIX.nwo"
MPLBACKEND=Agg python3 butane.py
```

Sustituye `butane.py` por el script del proyecto elegido. Ejecutarlo desde su carpeta es importante: varios scripts leen geometrías con rutas relativas. `MPLBACKEND=Agg` permite guardar gráficas sin requerir una sesión gráfica.

El backend de Quantum ESPRESSO usado en `carbon_dioxide/quantum_expresso` requiere el ejecutable `pw.x`, sus pseudopotenciales y los parámetros que espera `ase.calculators.espresso.Espresso`; instalar ASE por sí solo no instala esos componentes.

## Solución de problemas

- **`distrobox` o `podman` no se encuentra:** instala las herramientas en el host real, no solo dentro del terminal Flatpak de VS Code.
- **`podman info` falla:** configura Podman rootless siguiendo las instrucciones de tu distribución.
- **`nwchem` no se encuentra:** instala NWChem con la opción elegida y activa el entorno Conda si corresponde.
- **ASE no encuentra o no puede lanzar NWChem:** comprueba `command -v nwchem` y vuelve a exportar `ASE_NWCHEM_COMMAND` en la misma sesión.
- **No se importa ASE o una dependencia:** activa `ase-nwchem` o ejecuta el `python3` del sistema si instalaste con APT.
- **APT no encuentra los paquetes:** habilita `universe` y ejecuta `sudo apt update`.
- **No se encuentra un archivo XYZ:** entra primero a la carpeta del proyecto que contiene ese archivo.
- **El cálculo tarda o no converge:** inspecciona la salida de ASE/NWChem; el costo depende del sistema molecular y de los parámetros del cálculo.

Para la física, metodología, resultados y limitaciones de un ejemplo particular, consulta el README de ese proyecto.