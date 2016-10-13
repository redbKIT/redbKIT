#!/bin/bash

# INSTALL MUMPS

rm -rf MUMPS

mkdir -p MUMPS

cd MUMPS

wget http://mumps.enseeiht.fr/MUMPS_5.0.1.tar.gz

tar -xzf MUMPS_5.0.1.tar.gz

cd MUMPS_5.0.1

cp ../../Makefile_MUMPS.inc Makefile.inc

cp ../../make_MatlabMUMPS.inc MATLAB/make.inc

make clean

export rootPath=$(cd ../../; pwd)

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/scratch/Libraries/Test_MUMPS_Installation/OpenBlas/install/lib/

make all

cd examples

make all

make d

./dsimpletest < input_simpletest_real
./zsimpletest < input_simpletest_cmplx

cd ../MATLAB

export mexPath=$(which mex)

make clean

make

cp dmumps.m dmumpsmex.m* initmumps.m mumps_help.m ../../../../
