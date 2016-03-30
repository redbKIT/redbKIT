setPath;

this_path = pwd;

cd FEM_library/


%% FEM
cd Core
geotrasf_prj;
cd ../

%% ADR Models
cd Models/ADR/
% Generate C-Mex assemblers from matlab code and compile
ADR_assembly2D_prj;
ADR_assembly3D_prj;

% Compile C assembler
mex ADR_assembler2D_C.c
mex ADR_assembler2D_C_omp.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
mex ADR_assembler3D_C_omp.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
mex ADR_assembler_C_omp.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"

cd ../

%% ADR Models
cd Models/CSM/

% Compile C assembler
mex CSM_assembler_C_omp.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
mex CSM_assembler_ExtForces.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"

cd(this_path)