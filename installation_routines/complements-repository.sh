#!/bin/bash
#Scritp para la instalacion de complementos desde el repositorio
#Version 2.0
#Date: 05/04/2025
#Distribuciones Debian y derivadas.
#Arquitectura x86
echo "Verificar version de la distribucion"
sudo apt-get update
echo "Componentes basicos"
sudo apt-get -y install build-essential fftw3-dev gfortran
sudo apt-get -y install liblapack-dev fftw-dev
sudo apt-get -y install subversion libopenmpi-dev libopenmpi1.6 python-dev
echo "Calculadores"
sudo apt-get -y install nwchem
sudo apt-get -y install quantum-espresso
sudo apt-get -y install abinit
echo "Visualizadores"
sudo apt-get -y install pymol
sudo apt-get -y install openbabel
sudo apt-get -y install avogadro
