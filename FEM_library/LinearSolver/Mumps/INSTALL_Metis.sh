#!/bin/bash

# INSTALL METIS

rootPath=$PWD

rm -rf Metis

mkdir -p Metis

cd Metis

wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz

tar -xvzf metis-5.1.0.tar.gz

mkdir -p metis-5.1.0-install

cd metis-5.1.0

make config prefix=$rootPath/Metis/metis-5.1.0-install/

make

make install
