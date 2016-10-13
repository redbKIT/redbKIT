#!/bin/bash

#INSTALL SCOTCH

rootPath=$PWD

mkdir -p Scotch

cd Scotch

wget https://gforge.inria.fr/frs/download.php/file/34618/scotch_6.0.4.tar.gz

tar -xvzf scotch_6.0.4.tar.gz

cd scotch_6.0.4

cd src

touch Makefile.inc

echo "EXE		=" >> Makefile.inc
echo "LIB		= .a" >> Makefile.inc
echo "OBJ		= .o" >> Makefile.inc

echo "MAKE		= make" >> Makefile.inc
echo "AR		= ar" >> Makefile.inc
echo "ARFLAGS		= -ruv" >> Makefile.inc
echo "CAT		= cat" >> Makefile.inc
echo "CCS		= gcc -fPIC" >> Makefile.inc
echo "CCP		= mpicc" >> Makefile.inc
echo "CCD		= gcc -fPIC" >> Makefile.inc
echo "CFLAGS		= -O3 -DCOMMON_FILE_COMPRESS_GZ -DCOMMON_PTHREAD -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_RENAME -DSCOTCH_PTHREAD -Drestrict=__restrict -DIDXSIZE32 -DINTSIZE32" >> Makefile.inc
echo "CLIBFLAGS	=" >> Makefile.inc
echo "LDFLAGS		= -lz -lm -lrt -pthread" >> Makefile.inc
echo "CP		= cp" >> Makefile.inc
echo "LEX		= flex -Pscotchyy -olex.yy.c" >> Makefile.inc
echo "LN		= ln" >> Makefile.inc
echo "MKDIR		= mkdir" >> Makefile.inc
echo "MV		= mv" >> Makefile.inc
echo "RANLIB		= ranlib" >> Makefile.inc
echo "YACC		= bison -pscotchyy -y -b y" >> Makefile.inc

make

make esmumps 
