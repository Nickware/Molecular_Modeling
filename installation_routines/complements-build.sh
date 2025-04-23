#!/bin/bash
# Script de instalación de VMD y Quantum ESPRESSO
# Versión 1.3 – Fecha: 20/04/2025
# Para distribuciones Debian y derivadas – Arquitectura x86_64

set -e  # Detiene ejecución si ocurre algún error

# === CONFIGURACIÓN ===
VMD_VERSION="1.9.1"
INSTALL_DIR="/opt"
LOG_ESTADO="estado_instalacion.log"

# Limpiar log previo si existe
> "$LOG_ESTADO"

# === FUNCIONES ===

# Verificar si VMD se instaló correctamente
check_vmd() {
    echo -e "\n Verificando instalación de VMD..."
    if command -v vmd >/dev/null 2>&1; then
        echo " VMD está disponible en el sistema." | tee -a "$LOG_ESTADO"
    else
        echo " VMD no se encuentra en el PATH. Verifica la instalación o el .bashrc." | tee -a "$LOG_ESTADO"
    fi
}

# Verificar si Quantum ESPRESSO se instaló correctamente
check_qe() {
    echo -e "\n Verificando instalación de Quantum ESPRESSO..."
    if [ -x "$INSTALL_DIR/espresso/bin/pw.x" ]; then
        echo " Quantum ESPRESSO fue instalado exitosamente y 'pw.x' está presente." | tee -a "$LOG_ESTADO"
    else
        echo " No se encontró 'pw.x' en $INSTALL_DIR/espresso/bin. Revisa la compilación." | tee -a "$LOG_ESTADO"
    fi
}

# === INICIO DEL SCRIPT ===

echo " Iniciando instalación de herramientas científicas..."

# Actualizar paquetes
echo " Actualizando lista de paquetes..."
sudo apt-get update

# Instalar dependencias necesarias
DEPS=(build-essential gfortran wget tar make xorg libglu1-mesa-dev)
for pkg in "${DEPS[@]}"; do
    sudo apt-get -y install "$pkg"
done

# === INSTALACIÓN DE VMD ===
echo " Instalando VMD $VMD_VERSION..."

VMD_TAR="vmd-${VMD_VERSION}.bin.LINUXAMD64.opengl.tar.gz"
VMD_URL="http://www.ks.uiuc.edu/Research/vmd/vmd-${VMD_VERSION}/files/final/$VMD_TAR"

wget -c "$VMD_URL"
tar -zxvf "$VMD_TAR"
cd "vmd-${VMD_VERSION}"
./configure LINUXAMD64
cd src
sudo make install
cd ../..

echo " VMD instalado."
check_vmd

# === INSTALACIÓN DE QUANTUM ESPRESSO ===
echo "⚛ Instalación de Quantum ESPRESSO..."
echo " Pega la URL del archivo .tar.gz de Quantum ESPRESSO (GitLab):"
read -p "URL (ej: https://gitlab.com/...): " QE_URL

if [[ -z "$QE_URL" ]]; then
    echo " No se proporcionó una URL. Abortando..."
    exit 1
fi

QE_TAR=$(basename "$QE_URL")
QE_FOLDER="espresso"

wget -c "$QE_URL" -O "$QE_TAR"

if [[ $? -ne 0 ]]; then
    echo " Error al descargar el archivo. Verifica que la URL sea correcta."
    exit 1
fi

tar -xzf "$QE_TAR"
EXTRACTED_FOLDER=$(tar -tzf "$QE_TAR" | head -1 | cut -f1 -d"/")

cd "$EXTRACTED_FOLDER"
./configure
make all
cd ..
sudo mv "$EXTRACTED_FOLDER" "$INSTALL_DIR/$QE_FOLDER"

# Agregar al PATH del usuario
USER_HOME=$(eval echo ~${SUDO_USER})
echo "export PATH=\$PATH:$INSTALL_DIR/$QE_FOLDER/bin" >> "$USER_HOME/.bashrc"
source "$USER_HOME/.bashrc"

check_qe

echo -e "\n Estado de la instalación:"
cat "$LOG_ESTADO"

echo -e "\n Todo listo. Reinicia la terminal o ejecuta 'source ~/.bashrc' para usar los comandos instalados."
