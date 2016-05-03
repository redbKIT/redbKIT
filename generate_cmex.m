function generate_cmex(hasOpenMP)

if nargin < 1 || isempty(hasOpenMP)
    hasOpenMP = 0;
end

if hasOpenMP
    fprintf('\nCompiling with openmp enabled\n');
else
    fprintf('\nCompiling without openmp\n');
end

%setPath;

this_path = pwd;

cd FEM_library/

%% FEM
cd Models/

%% ADR Model
cd ADR/
% Generate C-Mex assemblers from matlab code and compile
%ADR_assembly2D_prj;
%ADR_assembly3D_prj;

% Compile C assembler
% to check whether openmp is available:
%   $ echo |cpp -fopenmp -dM |grep -i open
if hasOpenMP
    mex ADR_assembler_C_omp.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
else
    mex ADR_assembler_C_omp.c
end

if hasOpenMP
    mex Mass_assembler_C_omp.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
else
    mex Mass_assembler_C_omp.c
end

cd ../

%% CSM Model
cd CSM/

% Compile C assembler
if hasOpenMP
    mex CSM_assembler_C_omp.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
    mex CSM_assembler_ExtForces.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
    
else
    mex CSM_assembler_C_omp.c
    mex CSM_assembler_ExtForces.c
end

cd ../

%% CFD Model
cd CFD/

% Compile C assembler
if hasOpenMP
    mex CFD_assembler_C_omp.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
    mex CFD_assembler_ExtForces.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
    
else
    mex CFD_assembler_C_omp.c
    mex CFD_assembler_ExtForces.c
end

cd(this_path)

end