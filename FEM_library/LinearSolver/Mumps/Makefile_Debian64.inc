#
#  This file is part of MUMPS 5.0.1, released
#  on Thu Jul 23 17:08:29 UTC 2015
#
# These settings for a PC under Debian/linux with standard packages :
# metis (parmetis), scotch (ptscotch), openmpi, gfortran

# packages installation : 
# apt-get install libmetis-dev libparmetis-dev libscotch-dev libptscotch-dev libatlas-base-dev openmpi-bin libopenmpi-dev  

# Begin orderings
LSCOTCHDIR = /usr/lib/scotch-int32/lib
#ISCOTCH   = -I/usr/include/scotch # only needed for ptscotch

#LSCOTCH   = -L$(LSCOTCHDIR) -lptesmumps -lptscotch -lptscotcherr
LSCOTCH   = -L$(LSCOTCHDIR) -lesmumps -lscotch -lscotcherr

LPORDDIR = $(topdir)/PORD/lib/
IPORD    = -I$(topdir)/PORD/include/
LPORD    = -L$(LPORDDIR) -lpord

#LMETISDIR = /usr/lib 
#LMETISDIR = /usr/scratch/Dahlia/Code_repo/metis-5.0.2/libmetis
LMETISDIR = /usr/scratch/Libraries/Metis/metis-5.1.0_INT32-install/lib
#IMETIS    = -I/usr/include/parmetis
#IMETIS    = -I/usr/include/metis
IMETIS    = -I/usr/scratch/Libraries/Metis/metis-5.1.0_INT32-install/include
#IMETIS     = -I/usr/scratch/Dahlia/Code_repo/metis-5.0.2/include

# LMETIS    = -L$(LMETISDIR) -lparmetis
LMETIS    = -L$(LMETISDIR) -lmetis

# Corresponding variables reused later
#ORDERINGSF = -Dmetis -Dpord -Dparmetis -Dscotch -Dptscotch
ORDERINGSF = -Dmetis -Dpord -Dscotch
ORDERINGSC  = $(ORDERINGSF)
LORDERINGS = $(LMETIS) $(LPORD) $(LSCOTCH)
IORDERINGSF = $(ISCOTCH)
IORDERINGSC = $(IMETIS) $(IPORD) $(ISCOTCH)
# End orderings
################################################################################

PLAT    =
LIBEXT  = .a
OUTC    = -o 
OUTF    = -o
RM = /bin/rm -f
CC = gcc -fPIC
FC = gfortran -fPIC
FL = gfortran -fPIC
AR = ar vr 
RANLIB = ranlib



INCSEQ = -I$(topdir)/libseq
LIBSEQ  =  -L$(topdir)/libseq -lmpiseq

#LIBBLAS = -lblas
LIBBLAS = /usr/lib/openblas-base/libblas.so.3
# LIBOTHERS = -lpthread -lblas -lgfortran
LIBOTHERS = -lpthread /usr/lib/openblas-base/libblas.so.3 -lgfortran

#Preprocessor defs for calling Fortran from C (-DAdd_ or -DAdd__ or -DUPPER)
CDEFS   = -DAdd_

#Begin Optimized options
# OPTF    = -O3  -DALLOW_NON_INIT -fopenmp
OPTF    = -O3  -DALLOW_NON_INIT
# OPTL    = -O -fopenmp
OPTL    = -O
# OPTC    = -O3 -fopenmp
OPTC    = -O3
#End Optimized options
INCS = $(INCSEQ)
LIBS = $(LIBSEQ)
LIBSEQNEEDED = libseqneeded
