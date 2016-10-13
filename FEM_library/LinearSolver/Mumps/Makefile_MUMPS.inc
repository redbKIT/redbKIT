#
#  This file is part of MUMPS 5.0.1, released
#  on Thu Jul 23 17:08:29 UTC 2015
#
# These settings for a PC under Debian/linux with standard packages :
# metis (parmetis), scotch (ptscotch), openmpi, gfortran

# packages installation : 
# apt-get install libmetis-dev libparmetis-dev libscotch-dev libptscotch-dev libatlas-base-dev openmpi-bin libopenmpi-dev  

# Begin orderings
LSCOTCHDIR =  $(rootPath)/Scotch/scotch_6.0.4/lib
# LSCOTCHDIR = /usr/scratch/Libraries/Test_MUMPS_Installation/Scotch/scotch_6.0.4/lib

LSCOTCH   = -L$(LSCOTCHDIR) -lesmumps -lscotch -lscotcherr

LPORDDIR = $(topdir)/PORD/lib/
IPORD    = -I$(topdir)/PORD/include/
LPORD    = -L$(LPORDDIR) -lpord

LMETISDIR =  $(rootPath)/Metis/metis-5.1.0-install/lib
IMETIS    = -I$(rootPath)/Metis/metis-5.1.0-install/include

LMETIS    = -L$(LMETISDIR) -lmetis

# Corresponding variables reused later

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

LIBBLAS = -L$(rootPath)/OpenBlas/install/lib -lopenblas

LIBOTHERS = -lpthread

#Preprocessor defs for calling Fortran from C (-DAdd_ or -DAdd__ or -DUPPER)
CDEFS   = -DAdd_

#Begin Optimized options
OPTF    = -O3  -DALLOW_NON_INIT
OPTL    = -O
OPTC    = -O3
#End Optimized options

INCS = $(INCSEQ)
LIBS = $(LIBSEQ)
LIBSEQNEEDED = libseqneeded
