# It is possible to generate a MATLAB or an Octave interface thanks to
# the Octave MEX file compatibility. Comment/uncomment the lines below
# depending on whether you want to generate the MATLAB or the Octave
# interface

# To generate the MATLAB interface uncomment the following line
# ( the use of -largeArrayDims is necessary to work with sparse
# matrices since R2006b)
MEX         = /software/base/matlab/R2014a/bin/mex -g -largeArrayDims

# To generate the Octave interface uncomment the following line
# MEX         = mkoctfile -g --mex

# Main MUMPS_DIR
MUMPS_DIR = /srv/share/malc1/home/negri/Libraries/MUMPS/Sources

# Orderings (see main Makefile.inc file from MUMPS)
LMETISDIR = /home-malc1/negri/Libraries/Metis/metis-5.1.0_INT32-install/lib
LMETIS     = -L$(LMETISDIR) -lmetis
LPORDDIR   = $(MUMPS_DIR)/PORD/lib
LPORD      = -L$(LPORDDIR) -lpord
# LORDERINGS = $(LMETIS) $(LPORD)


LSCOTCHDIR = /home-malc1/negri/Libraries/Scotch/scotch_6.0.4/lib
LSCOTCH   = -L$(LSCOTCHDIR) -lesmumps -lscotch -lscotcherr

LORDERINGS = $(LMETIS) $(LPORD) $(LSCOTCH)

# Fortran runtime library
# Please find out the path and name of your
# Fortran runtime, examples below:
# g95:
# LIBFORT = /usr/lib/libf95.a /usr/lib/libgcc.a 
# Intel:
# LIBFORT = /opt/intel80/lib/libifcore.a /opt/intel80/lib/libifport.a /opt/intel80/lib/libirc.a
# PGI:
# LIBFORT = -L/usr/local/pgi/linux86/5.2/lib -llapack -lblas -lpgf90 -lpgc -lpgf90rtl -lpgftnrtl -lpgf902 -lpgf90_rpm1 -lpghpf2
# SGI 32-bit
# LIBFORT = -L/usr/lib32 -lblas -L/usr/lib32/mips4 -lfortran
# Sun
# LIBFORT = -L/opt2/SUNWspro7/lib -lsunperf -lfminvai -lfai2 -lfsu -lfmaxvai -lfmaxlai -lfai -lfsumai -lLIBFORT = /usr/local/lib/libgfortran.a

# gfortran
  LIBFORT = -lgfortran

# BLAS library:
# LIBBLAS = -L/usr/lib/atlas -lblas
# LIBBLAS = -lsunperf -lf77compat
# LIBBLAS = -lblas
# LIBBLAS = /usr/lib/openblas-base/libblas.so
# LIBBLAS = -I/usr/scratch/Libraries/OpenBlas/install/include -L/usr/scratch/Libraries/OpenBlas/install/lib/ -lmyopenblas
LIBBLAS = -I/home-malc1/negri/Libraries/OpenBlas/install/include -L/home-malc1/negri/Libraries/OpenBlas/install/lib -lmyopenblas


# extra options passed via mex command
# Add -DINTSIZE64 if MUMPS was compiled with 64-bit integers (BLAS
# library should then have 64-bit integers as well)
#OPTC    = -O CFLAGS="-fopenmp -fPIC" LDFLAGS="-fopenmp -fPIC"
OPTC    = -O
