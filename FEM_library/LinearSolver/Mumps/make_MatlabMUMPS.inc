# It is possible to generate a MATLAB or an Octave interface thanks to
# the Octave MEX file compatibility. Comment/uncomment the lines below
# depending on whether you want to generate the MATLAB or the Octave
# interface

# To generate the MATLAB interface uncomment the following line
# ( the use of -largeArrayDims is necessary to work with sparse
# matrices since R2006b)
MEX         = $(mexPath) -g -largeArrayDims

# To generate the Octave interface uncomment the following line
# MEX         = mkoctfile -g --mex

# Main MUMPS_DIR
MUMPS_DIR = $(rootPath)/MUMPS/MUMPS_5.0.1

# Orderings (see main Makefile.inc file from MUMPS)
LMETISDIR = $(rootPath)/Metis/metis-5.1.0-install/lib
LMETIS     = -L$(LMETISDIR) -lmetis
LPORDDIR   = $(MUMPS_DIR)/PORD/lib
LPORD      = -L$(LPORDDIR) -lpord
# LORDERINGS = $(LMETIS) $(LPORD)


LSCOTCHDIR = $(rootPath)/Scotch/scotch_6.0.4/lib
LSCOTCH   = -L$(LSCOTCHDIR) -lesmumps -lscotch -lscotcherr

LORDERINGS = $(LMETIS) $(LPORD) $(LSCOTCH)

# gfortran
LIBFORT = -lgfortran

# BLAS library:
#LIBBLAS = -I$(rootPath)/OpenBlas/install/include -L$(rootPath)/OpenBlas/install/lib -lopenblas
LIBBLAS = -I$(rootPath)/OpenBlas/install/include -L$(rootPath)/OpenBlas/install/lib -lmyopenblas

# extra options passed via mex command
# Add -DINTSIZE64 if MUMPS was compiled with 64-bit integers (BLAS
# library should then have 64-bit integers as well)
#OPTC    = -O CFLAGS="-fopenmp -fPIC" LDFLAGS="-fopenmp -fPIC"
OPTC    = -O
