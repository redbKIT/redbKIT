# It is possible to generate a MATLAB or an Octave interface thanks to
# the Octave MEX file compatibility. Comment/uncomment the lines below
# depending on whether you want to generate the MATLAB or the Octave
# interface

# To generate the MATLAB interface uncomment the following line
# ( the use of -largeArrayDims is necessary to work with sparse
# matrices since R2006b)
MEX         = /usr/local/MATLAB/R2016a/bin/mex -g -largeArrayDims

# To generate the Octave interface uncomment the following line
# MEX         = mkoctfile -g --mex

# Main MUMPS_DIR
MUMPS_DIR = /usr/scratch/Libraries/MUMPS/SEQ/Sources

# Orderings (see main Makefile.inc file from MUMPS)
LMETISDIR = /usr/scratch/Libraries/Metis/metis-5.1.0_INT32-install/lib
# LMETISDIR = /usr/scratch/Dahlia/Code_repo/metis-5.0.2/libmetis
LMETIS     = -L$(LMETISDIR) -lmetis
LPORDDIR   = $(MUMPS_DIR)/PORD/lib
LPORD      = -L$(LPORDDIR) -lpord
# LORDERINGS = $(LMETIS) $(LPORD)

#LSCOTCHDIR = /usr/lib/scotch_int32
LSCOTCHDIR =/usr/lib/scotch-int32/lib
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
  LIBFORT = /usr/lib/gcc/x86_64-linux-gnu/4.9.2/libgfortran.so

# BLAS library:
# LIBBLAS = -L/usr/lib/atlas -lblas
# LIBBLAS = -lsunperf -lf77compat
# LIBBLAS = -lblas
LIBBLAS = /usr/lib/openblas-base/libblas.so
#LIBBLAS = /home/jylexcel/libs_courge/libgoto_opteronp-r1.26.a
# LIBBLAS = /home/jylexcel/libs_courge/libblas.a  # -fPIC missing
# LIBBLAS = /usr/lib/libblas.so

# extra options passed via mex command
# Add -DINTSIZE64 if MUMPS was compiled with 64-bit integers (BLAS
# library should then have 64-bit integers as well)
#OPTC    = -O CFLAGS="-fopenmp -fPIC" LDFLAGS="-fopenmp -fPIC"
OPTC    = -O
