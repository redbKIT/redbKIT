function [ R ] = CSM_overlapping_DD( MESH, n_subdom, overlap )
%CSM_OVERLAPPING_DD builds subdomains with overlap 
%
%   INPUT: dim: 2 for 2D, 3 for 3D
%        n_subdom: number of subdomains
%        overlap: number of overlap layers, 1 is suggested

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

vertices = MESH.vertices;
elements = MESH.elements(1:MESH.dim+1,:);% P1 elements
dim      = MESH.dim;
I        = MESH.internal_dof_c;
I_all    = MESH.internal_dof;
nln      = MESH.numElemDof;
nodes    = MESH.nodes;
elements_fem = MESH.elements(1:nln,:);% P1 elements

[subdom] = geometric_domain_decomposition(vertices, elements, dim, n_subdom, overlap, 1, 'Figures', elements_fem);

%% restrict subdomains to internal vertices for each displacement component
n_component = dim;
for k = 1 : n_component
      for i = 1 : n_subdom
            subdom_I{i,k} = intersect(subdom{i}, I{k});
            if size(subdom_I{i,k},1) > size(subdom_I{i,k},2)
                  subdom_I{i,k} = subdom_I{i,k}';
            end
      end
end

%% build subdomains restriction/prolongation operators
nov     = size(nodes,2);
nov_tot = dim*nov;

parfor i = 1 : n_subdom
   
    subdom_I_all = [];
 
    for k = 1 : n_component
        subdom_I_all = [subdom_I_all nov*(k-1)+subdom_I{i,k}];
    end
 
    tmp =  zeros(nov_tot,1);
    tmp(subdom_I_all) = 1;
    tmp2  =  tmp(I_all);    
    R{i}  = find(tmp2);
    
end

check_n_subd = length(R);
fprintf('\n%d subdomains and restriction/prolongation operators built ---\n',check_n_subd);
        
return