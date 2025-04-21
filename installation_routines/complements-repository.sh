#!/bin/bash
# Script para la instalación de complementos desde el repositorio
# Versión 2.1 - 20/04/2025
# Para distribuciones Debian y derivadas, arquitectura x86_64

# Limpiar logs anteriores
rm -f paquetes_instalados.log paquetes_fallidos.log

# Verificar distribución
. /etc/os-release
echo " Distribución detectada: $PRETTY_NAME"

# Verificar arquitectura
arch=$(uname -m)
if [[ "$arch" != "x86_64" ]]; then
    echo " Advertencia: Arquitectura detectada: $arch. El script fue diseñado para x86_64."
fi

# Actualizar repositorios
echo " Actualizando lista de paquetes..."
sudo apt-get update || { echo " Falló apt-get update"; exit 1; }

# Función para verificar e instalar paquetes
install_pkg() {
    pkg="$1"
    echo -e "\n Verificando disponibilidad de $pkg..."

    if apt-cache show "$pkg" >/dev/null 2>&1; then
        echo " $pkg está disponible. Instalando..."
        if sudo apt-get -y install "$pkg"; then
            echo " $pkg instalado correctamente."
            echo "$pkg (instalado correctamente)" >> paquetes_instalados.log
        else
            echo " Error durante la instalación de $pkg" >&2
            echo "$pkg (error durante instalación)" >> paquetes_fallidos.log
        fi
    else
        echo " $pkg no está disponible en los repositorios." >&2
        echo "$pkg (no disponible)" >> paquetes_fallidos.log
        echo " Sugerencias de paquetes similares:"
        apt-cache search "$pkg" | head -5
    fi
}

# Componentes básicos
echo -e "\n Instalando componentes básicos..."
for pkg in build-essential fftw3-dev gfortran liblapack-dev fftw-dev; do
    install_pkg "$pkg"
done

# Calculadores científicos
echo -e "\n Instalando calculadores científicos..."
for pkg in nwchem quantum-espresso abinit; do
    install_pkg "$pkg"
done

# Visualizadores moleculares
echo -e "\n Instalando visualizadores moleculares..."
for pkg in pymol openbabel avogadro; do
    install_pkg "$pkg"
done

# Mensaje final
echo -e "\n Instalación finalizada."
echo " Revisa 'paquetes_instalados.log' y 'paquetes_fallidos.log' para más detalles."

