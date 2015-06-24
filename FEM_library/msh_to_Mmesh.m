function [vertices, boundaries, elements] = msh_to_Mmesh(filename, dimension)
%GMSH_TO_MMESH import a msh mesh file .msh and transform it to a 'pdetool'
%like Matlab format
%
%   [VERTICES, BOUNDARIES, ELEMENTS] = GMSH_TO_MMESH(FILENAME, DIMENSION)

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

fprintf('\nLoading mesh from file %s.msh ... ',filename);

if nargin < 2
      dimension = 2;
end

time_load = tic;

if dimension == 2
      mesh_filename    =  strcat(filename,'.msh');
      mesh             =  read_msh(mesh_filename);
      vertices         =  mesh.NODES(1:2,:);
      elements         =  mesh.ELEMENTS{2};
      tmp_boundaries   =  mesh.ELEMENTS{1};
      
      boundaries(1:2,:)  = tmp_boundaries(1:2,:);
      boundaries(3:4,:)  = 0*tmp_boundaries(1:2,:);
      boundaries(5,:)    = tmp_boundaries(3,:);
      boundaries(6,:)    = ones(size(tmp_boundaries(3,:)));
      boundaries(7,:)    = 0*tmp_boundaries(1,:);
      
elseif dimension == 3  
      mesh_filename    =  strcat(filename,'.msh');
      mesh             =  read_msh(mesh_filename);
      vertices         =  mesh.NODES;
      elements         =  mesh.ELEMENTS{4};
      tmp_boundaries   =  mesh.ELEMENTS{2};
      
      boundaries([1 2 3 12],:)   = tmp_boundaries([1 2 3 4],:);
end

time_load = toc(time_load);

fprintf('done in %3.2f s \n\n',time_load);


return