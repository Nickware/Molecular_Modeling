# Instalación y ejecución del ejemplo de butano con Distrobox

Esta guía crea un contenedor Ubuntu 24.04 y prepara en él NWChem, ASE y las bibliotecas Python que usa `butane.py`. Separa los comandos del sistema anfitrión de los que se ejecutan dentro del contenedor.

## 1. Comprobar el sistema anfitrión

Distrobox necesita Linux y un motor de contenedores. Se recomienda Podman en modo rootless; Docker también es compatible si está instalado y configurado para que el usuario pueda ejecutarlo sin `sudo`.

Ejecuta en una terminal del sistema anfitrión:

```bash
cat /etc/os-release
```

Si el resultado identifica un runtime como `Freedesktop SDK` en lugar de tu distribución (por ejemplo, al usar VS Code instalado como Flatpak), esa terminal está dentro del sandbox de la aplicación. Instala y administra Distrobox desde una terminal del anfitrión, fuera de ese sandbox. El terminal de VS Code no necesariamente ve sus ejecutables aunque estén instalados en el host.

Instala Distrobox y Podman usando el gestor de paquetes de tu distribución. Ejemplos para sistemas que incluyen ambos paquetes en sus repositorios:

```bash
# Ubuntu 22.10 o posterior
sudo apt update
sudo apt install -y distrobox podman

# Fedora
sudo dnf install -y distrobox podman

# Arch Linux
sudo pacman -Syu distrobox podman
```

En otras versiones o distribuciones, sigue la [guía oficial de instalación y compatibilidad de Distrobox](https://distrobox.it/installation/). No ejecutes a ciegas un instalador remoto con privilegios de administrador.

Comprueba que Distrobox y Podman están instalados y que el motor funciona como usuario normal:

```bash
distrobox --version
podman --version
podman info
```

`podman info` debe finalizar sin errores. Si eliges Docker, comprueba `docker info` y configura los permisos según la [documentación oficial](https://docs.docker.com/engine/install/linux-postinstall/). Distrobox es un entorno integrado con el host, no un sandbox de seguridad: el contenedor puede acceder a tu directorio personal.

## 2. Crear y entrar al contenedor

El repositorio debe estar en tu directorio personal, que Distrobox monta en el contenedor. Si todavía no lo tienes, clónalo en el host antes de crear el contenedor:

```bash
mkdir -p ~/gitHub
cd ~/gitHub
git clone https://github.com/Nickware/Molecular_Modeling.git MolecularModeling
```

Si ya tienes el proyecto en otro directorio dentro de `$HOME`, no necesitas volver a clonarlo.

En el host, crea el contenedor y entra:

```bash
distrobox create --name butane-dft --image ubuntu:24.04
distrobox enter --name butane-dft
```

La primera creación descarga la imagen y puede tardar. Los pasos restantes se ejecutan dentro de la sesión Ubuntu del contenedor.

## 3. Instalar ASE y NWChem

Dentro del contenedor, comprueba que tienes permisos para instalar paquetes:

```bash
whoami
sudo -v
```

Escoge **una** de estas dos alternativas. No mezcles los paquetes Python de APT con los de Conda: cada opción mantiene ASE y sus dependencias en un solo entorno.

### Opción A: Conda-forge (entorno aislado)

Instala Miniforge dentro del contenedor. El siguiente comando elige el instalador según la arquitectura Linux:

```bash
sudo apt update
sudo apt install -y ca-certificates curl
curl -L "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$(uname -m).sh" -o /tmp/Miniforge3.sh
bash /tmp/Miniforge3.sh -b -p "$HOME/miniforge3"
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda create -y --override-channels -c conda-forge -n butane-dft python=3.12 ase nwchem numpy matplotlib
conda activate butane-dft
```

Miniforge es el instalador de conda-forge. En sesiones futuras, carga `conda.sh` y activa `butane-dft` antes de usar Python.

### Opción B: repositorios de Ubuntu 24.04

Ubuntu publica `nwchem` y `python3-ase` en el componente `universe`. Instala también las dependencias Python usadas por el script:

```bash
sudo apt update
sudo apt install -y nwchem python3-ase python3-numpy python3-matplotlib
```

Comprueba que APT encontró ambos paquetes:

```bash
apt-cache policy nwchem python3-ase
```

Si no hay versión candidata, habilita `universe` y vuelve a instalar:

```bash
sudo apt install -y software-properties-common
sudo add-apt-repository -y universe
sudo apt update
sudo apt install -y nwchem python3-ase python3-numpy python3-matplotlib
```

Con esta opción ejecuta el Python del sistema (`python3`), no un entorno virtual que no vea los módulos instalados por APT.

Si `sudo` no está instalado o el usuario no tiene privilegios, no eleves permisos ejecutando `sudo distrobox`. Sal de la sesión y consulta `distrobox enter --help` para el modo root soportado por el motor configurado, o corrige la configuración de la imagen según la documentación de Distrobox.

## 4. Configurar el calculador NWChem de ASE

NWChem es el ejecutable que realiza el cálculo; ASE es la interfaz Python que genera sus archivos de entrada, lo lanza y lee sus resultados. Comprueba que el ejecutable esté en el `PATH` del entorno elegido:

```bash
command -v nwchem
```

Configura el comando que ASE usará. ASE sustituye `PREFIX` por el nombre base de sus archivos:

```bash
export ASE_NWCHEM_COMMAND="nwchem PREFIX.nwi > PREFIX.nwo"
```

La clase usada por el proyecto es `ase.calculators.nwchem.NWChem`. El comando configurado permite que `NWChem(xc='PBE')` en `butane.py` encuentre y ejecute NWChem. Comprueba que ASE puede importar el calculador:

```bash
python - <<'PY'
from ase import __version__
from ase.calculators.nwchem import NWChem

print("ASE:", __version__)
print("Calculador:", NWChem(xc="PBE", basis="3-21G"))
PY
```

Si abres otra sesión del contenedor, vuelve a activar el entorno Conda (si elegiste esa opción) y a exportar `ASE_NWCHEM_COMMAND` antes de ejecutar el script.

## 5. Prueba de instalación: cálculo pequeño con H2

Esta prueba valida Python, ASE, la configuración del comando y el ejecutable NWChem en un cálculo real. Los archivos se guardan en un directorio temporal:

```bash
python - <<'PY'
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

La prueba debe terminar sin traceback e imprimir una energía numérica. Si falla, revisa primero `command -v nwchem`, que `ASE_NWCHEM_COMMAND` esté exportada y que el entorno Python elegido esté activo. Conserva el mensaje completo de NWChem si el problema continúa.

## 6. Ejecutar el cálculo de butano

Desde el contenedor, cambia a la carpeta del ejemplo. Ajusta la ruta si el repositorio está en otro lugar dentro de tu directorio personal:

```bash
cd ~/gitHub/MolecularModeling/butane

# Solo para la opción Conda:
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate butane-dft

export ASE_NWCHEM_COMMAND="nwchem PREFIX.nwi > PREFIX.nwo"
MPLBACKEND=Agg python3 butane.py
```

El script lee `butane.xyz`, optimiza la geometría con ASE/BFGS y NWChem, calcula las frecuencias vibracionales e imprime la energía, la fuerza residual y el diagnóstico de modos imaginarios. `MPLBACKEND=Agg` permite guardar la gráfica aunque la sesión no tenga interfaz gráfica.

Al terminar, los resultados principales están en la carpeta `butane` del repositorio:

- `butane_opt.xyz`: geometría optimizada.
- `optimization.traj`: trayectoria BFGS.
- `diagnostico_frecuencias.png`: gráfica de frecuencias.
- La terminal: energía, fuerza máxima residual y número de modos imaginarios.

ASE también puede crear archivos de caché vibracional. No borres los resultados que quieras conservar: el script auxiliar `butane.sh` elimina varios archivos de salida.

## 7. Volver a entrar y detener el contenedor

Para otra sesión, ejecuta desde el host:

```bash
distrobox enter --name butane-dft
```

Para detener el contenedor sin borrarlo:

```bash
distrobox stop butane-dft
```

## Solución de problemas

- **`distrobox` o `podman` no se encuentra:** confirma que lo instalaste en el host real, no solamente dentro del terminal Flatpak de VS Code.
- **`nwchem` no se encuentra:** confirma que instalaste NWChem en el contenedor y que activaste `butane-dft` si usaste Conda.
- **ASE falla al lanzar NWChem:** comprueba que `ASE_NWCHEM_COMMAND` está exportada con el valor indicado y que `command -v nwchem` funciona en esa misma sesión.
- **Falta un módulo Python:** activa el entorno Conda, o confirma que ejecutas `python3` del sistema si instalaste ASE por APT.
- **`podman info` falla:** configura Podman rootless según la documentación de tu distribución. Distrobox requiere un motor funcional.
- **`sudo` no se encuentra en el contenedor:** la imagen o su inicialización no proporciona la configuración esperada; consulta `distrobox enter --help` y la guía oficial antes de instalar paquetes como usuario sin privilegios.
- **APT no encuentra `nwchem`:** confirma `apt-cache policy nwchem`; en Ubuntu habilita `universe` con los comandos de arriba.
- **La prueba de H2 falla:** no ejecutes todavía el cálculo completo. Resuelve primero el error de NWChem/ASE indicado por la salida.
- **No se encuentra `butane.xyz`:** ejecuta el script desde `butane/`, porque el código usa una ruta relativa al directorio actual.
- **La optimización tarda o no converge:** revisa la salida de ASE/NWChem. El script deja algunos parámetros electrónicos en sus valores predeterminados y el costo depende de los recursos del equipo.

Para la física, la metodología y las limitaciones científicas del cálculo, consulta el [README del proyecto de butano](../Readme.md).