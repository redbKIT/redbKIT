setPath;

cd FEM_library/

% Generate C-Mex assemblers from matlab code and compile
ADR_assembly2D_prj;
ADR_assembly3D_prj;

% Compile C assembler
mex ADR_assembler2D_C.c
mex ADR_assembler2D_C_omp.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
mex ADR_assembler3D_C_omp.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"

cd ../