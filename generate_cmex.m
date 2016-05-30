function generate_cmex(hasOpenMP)
%GENERATE_CMEX compiles finite element assemblers written in C
%
%   GENERATE_CMEX(HASOPENMP) if HASOPENMP = 1, compiles with flag -fopenmp.
%   Default is HASOPENMP = 0. 
%
%   On Linux: to check whether openmp is available on your system, open a 
%   terminal and type:
%   $ echo |cpp -fopenmp -dM |grep -i open

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

if nargin < 1 || isempty(hasOpenMP)
    hasOpenMP = 0;
end


%% Compile C code

if hasOpenMP
    fprintf('\nCompiling with openmp enabled\n');
    Flags = 'CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"';

else
    fprintf('\nCompiling without openmp\n');
    Flags = '';
end

source_files{1} = {'FEM_library/Models/ADR/','ADR_assembler_C_omp.c'};
source_files{2} = {'FEM_library/Models/ADR/','Mass_assembler_C_omp.c'};
source_files{3} = {'FEM_library/Models/CSM/','CSM_assembler_C_omp.c'};
source_files{4} = {'FEM_library/Models/CSM/','CSM_assembler_ExtForces.c'};
source_files{5} = {'FEM_library/Models/CFD/','CFD_assembler_C_omp.c'};
source_files{6} = {'FEM_library/Models/CFD/','CFD_assembler_ExtForces.c'};
source_files{7} = {'FEM_library/Models/CSM/','CSM_ComputeStress_C_omp.c'};

for i = 1 : length(source_files)
    
    file_path = [pwd, '/', source_files{i}{1}];
    file_name = source_files{i}{2};
    
    mex_command = sprintf( 'mex %s%s %s -outdir %s', file_path, file_name, Flags, file_path);
    eval( mex_command );
    
end

%% Mexify some Matlab codes

here = pwd;
source_files{1} = {'FEM_library/Models/CSM/','mexify_CSM_Assembly_M'};

for i = 1 : length(source_files)
    
    file_path = [pwd, '/', source_files{i}{1}];
    file_name = source_files{i}{2};
        
    eval(['cd ', file_path]);
    eval( file_name );
    eval(['cd ', here]);
    
end

end