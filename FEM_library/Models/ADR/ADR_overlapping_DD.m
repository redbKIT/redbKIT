function [ R ] = ADR_overlapping_DD( MESH, n_subdom, overlap, n_aggregates )
%ADR_OVERLAPPING_DD builds restriction operators associated to mesh decompositions
%for CSM problems
%
%   [ R ] = ADR_OVERLAPPING_DD( MESH, N_SUBDOM, OVERLAP )
%   given a MESH data structure (see also buildMESH.m), the number of
%   subdomains N_SUBDOM and the overlap level OVERLAP (>= 1), returns a cell
%   array R of length N_SUBDOM. Each element of R (say R{i}) contains a
%   vector with the indices of the DOFs pertaining to the i-th subdomain.
%
%   see also geometric_domain_decomposition, metis_to_matlab,
%   AS_Preconditioner

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

if nargin < 4 || isempty(n_aggregates)
    compute_coarse_aggregates = false;
else
    compute_coarse_aggregates = true;
end

vertices     = MESH.vertices;
elements     = MESH.elements(1:MESH.dim+1,:);% P1 elements
dim          = MESH.dim;
I            = MESH.internal_dof;
nln          = MESH.numElemDof;
nodes        = MESH.nodes;
elements_fem = MESH.elements(1:nln,:);% P1 elements

[subdom, ~, A, A_elemC] = geometric_domain_decomposition(vertices, elements, dim, n_subdom, overlap, 1, 'Figures', elements_fem);

%% restrict subdomains to internal vertices
for i = 1 : n_subdom
    subdom_I{i} = intersect(subdom{i}, I);
    if size(subdom_I{i},1) > size(subdom_I{i},2)
        subdom_I{i} = subdom_I{i}';
    end
end


%% build subdomains restriction/prolongation operators
nov     = size(nodes,2);

parfor i = 1 : n_subdom
    tmp =  zeros(nov,1);
    tmp(subdom_I{i}) = 1;
    tmp2  =  tmp(I);    
    R{i}  = find(tmp2);
end

check_n_subd = length(R);
fprintf('\n%d subdomains and restriction/prolongation operators built ---\n',check_n_subd);
        
%% build coarse aggregation restriction/prolongation operator

if compute_coarse_aggregates
    
    [subdom_Coarse] = geometric_aggregates(A_elemC, vertices, elements, dim, n_aggregates, elements_fem(1:MESH.numElemDof,:));
    
    R{n_subdom+1}  = sparse(n_aggregates, nov);
    
    for i = 1 : n_aggregates
        R{n_subdom+1}(i, subdom_Coarse{i}) = 1;
    end
    
    R{n_subdom+1} = R{n_subdom+1}(:,I);
    
    fprintf('\n%d Coarse aggregates computed ---\n', n_aggregates);
    
end


return