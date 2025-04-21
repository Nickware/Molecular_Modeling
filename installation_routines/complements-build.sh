#!/bin/bash
#============================================================
# Script de instalación automática de herramientas científicas VMD y Quantum Espresso
# Versión: 1.0
# Fecha: 05/04/2025
# Compatible con: Distribuciones basadas en Debian (x86_64)
# Herramientas: VMD 1.9.1 y Quantum ESPRESSO 5.4.0
#============================================================

# ------------------------------
# Parámetros configurables
# ------------------------------
VMD_VERSION="1.9.1"
QE_VERSION="5.4.0"
INSTALL_DIR="/opt"
USER_HOME="/home/$SUDO_USER"

# ------------------------------
# Verificación del usurio root
# ------------------------------
if [[ $EUID -ne 0 ]]; then
  echo "Este script debe ejecutarse como root. Usa: sudo $0"
  exit 1
fi

# ------------------------------
# ACTUALIZACIÓN DE REPOSITORIOS
# ------------------------------
echo " Actualizando repositorios de paquetes..."
apt-get update -y

# ------------------------------
# Instalación de dependecias básicas
# ------------------------------
echo " Instalando dependencias requeridas para compilar..."
apt-get install -y build-essential gfortran wget tar \
  libfftw3-dev libblas-dev liblapack-dev libxc-dev \
  libopenmpi-dev openmpi-bin xauth libglu1-mesa

# ------------------------------
# Instalación VMD (decargar, construir e instalar)
# ------------------------------
echo " Descargando VMD $VMD_VERSION..."
cd /tmp
wget -c http://www.ks.uiuc.edu/Research/vmd/vmd-$VMD_VERSION/files/final/vmd-$VMD_VERSION.bin.LINUXAMD64.opengl.tar.gz

echo " Extrayendo VMD..."
tar -xzf vmd-$VMD_VERSION.bin.LINUXAMD64.opengl.tar.gz
cd vmd-$VMD_VERSION
./configure LINUXAMD64
cd src
make install
cd /tmp

# ------------------------------
# Instalación Quantum Espresso (descargar, construir e instalar)
# ------------------------------
echo " Descargando Quantum ESPRESSO $QE_VERSION..."
wget -c http://qe-forge.org/gf/download/frsrelease/211/968/espresso-$QE_VERSION.tar.gz

echo " Extrayendo Quantum ESPRESSO..."
tar -xvf espresso-$QE_VERSION.tar.gz
mv espresso-$QE_VERSION $INSTALL_DIR/espresso

echo " Compilando Quantum ESPRESSO..."
cd $INSTALL_DIR/espresso
./configure
make all

# ------------------------------
# Configuración del entorno
# ------------------------------
echo " Configurando PATH para Quantum ESPRESSO..."
QE_PATH_LINE="export PATH=\$PATH:$INSTALL_DIR/espresso/bin"
if ! grep -q "$QE_PATH_LINE" "$USER_HOME/.bashrc"; then
  echo "$QE_PATH_LINE" >> "$USER_HOME/.bashrc"
  echo " PATH de Quantum ESPRESSO agregado a ~/.bashrc"
fi

# ------------------------------
# Finalización
# ------------------------------
echo " Instalación completada exitosamente."
echo " Reinicia tu terminal o ejecuta: source ~/.bashrc"

