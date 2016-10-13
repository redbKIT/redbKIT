#!/bin/bash

# INSTALL OpenBlas

rootPath=$PWD

mkdir -p OpenBlas

cd OpenBlas

wget http://github.com/xianyi/OpenBLAS/archive/v0.2.19.tar.gz

tar -xvzf v0.2.19.tar.gz

mkdir -p install

cd OpenBLAS-0.2.19

make FC=gfortran
make PREFIX=$rootPath/OpenBlas/install/ install 

cd ../install/lib

ln -s libopenblas.a libmyopenblas.a
