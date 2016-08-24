function make(hasOpenMP, sources)
%MAKE compiles finite element assemblers written in C
%
%   MAKE(HASOPENMP) if HASOPENMP = 1, compiles with flag -fopenmp.
%   Default is HASOPENMP = 0. 
%
%   MAKE(HASOPENMP, SOURCES) SOURCES is a vector of indices specifying the
%   source files to be compiled, e.g. MAKE(1, 8) compiles source Nr 8 with
%   OpenMP enabled.
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
dependencies{1} = {};
source_files{2} = {'FEM_library/Models/ADR/','Mass_assembler_C_omp.c'};
dependencies{2} = {};
source_files{3} = {'FEM_library/Models/CSM/','CSM_assembler_ExtForces.c'};
dependencies{3} = {};
source_files{4} = {'FEM_library/Models/CFD/','CFD_assembler_C_omp.c'};
dependencies{4} = {'../../Core/Tools.c'};
source_files{5} = {'FEM_library/Models/CFD/','CFD_assembler_ExtForces.c'};
dependencies{5} = {};
source_files{6} = {'FEM_library/Models/CSM/','CSM_assembler_C_omp.c'};
dependencies{6} = {'../../Core/Tools.c', 'MaterialModels/NeoHookeanMaterial.c',...
                   'MaterialModels/LinearElasticMaterial.c', 'MaterialModels/SEMMTMaterial.c', ...
                   'MaterialModels/StVenantKirchhoffMaterial.c',...
                   'MaterialModels/RaghavanVorpMaterial.c'};
source_files{7} = {'RB_library/Tools/RBF_interpolation/','RBF_evaluate_Fast.c'};
dependencies{7} = {};
source_files{8} = {'FEM_library/Models/ADR/','ADR_SUPGassembler_C_omp.c'};
dependencies{8} = {'../../Core/Tools.c'};

%Mexify = 0;               
if nargin < 2 || isempty( sources )
    sources = 1:length(source_files);
    %Mexify  = 1;
end
    
n_sources = length(sources); 

k = 0;

for i = sources
    
    k = k + 1;
    
    fprintf('\n ** Compiling source Nr. %d; %d of %d \n', i, k, n_sources)
    file_path = [pwd, '/', source_files{i}{1}];
    file_name = source_files{i}{2};
    
    all_dep = '';
    if ~isempty( dependencies{i} )
        for j = 1 : length( dependencies{i})
            
            this_dep =  sprintf('%s%s', file_path, dependencies{i}{j});
            all_dep = [all_dep, ' ', this_dep]; 
        end
        
    end
    
    mex_command = sprintf( 'mex %s%s %s %s -outdir %s', file_path, file_name, all_dep, Flags, file_path);
    eval( mex_command );
end

%% Mexify some Matlab codes

% here = pwd;
% %source_files_M{1} = {'FEM_library/Models/CSM/','mexify_CSM_Assembly_M'};
% %source_files_M{2} = {'FEM_library/Models/CFD/','mexify_SUPG'};
% 
% if Mexify
%     for i = 1 : length(source_files_M)
%         
%         file_path = [pwd, '/', source_files_M{i}{1}];
%         file_name = source_files_M{i}{2};
%         
%         eval(['cd ', file_path]);
%         eval( file_name );
%         rmdir('codegen', 's');
%         eval(['cd ', here]);
%         
%     end
% end

end