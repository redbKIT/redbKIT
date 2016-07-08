function [subdom, subdom_noOverlap, A, A_elemC] = geometric_domain_decomposition(vertices, elements, dim, n_subdom, overlap, visual, folder, elements_fem)
%GEOMETRIC_DOMAIN_DECOMPOSITION builds mesh subdomains with and without overlap 
%using Metis Library
%
%   [subdom, subdom_noOverlap] = GEOMETRIC_DOMAIN_DECOMPOSITION(vertices, elements, dim, n_subdom, overlap, visual, folder, elements_fem)
%
%   see also metis_to_matlab.m

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

if nargin < 6
      visual = 1;
end

if (visual && nargin < 7) || isempty(folder)
    folder = 'Figures_DD';
end

if nargin < 8 || isempty(elements_fem)
    elements_fem = elements;
    nln          = dim + 1;% for P1 FEM in dim=2,3
    elements_fem = elements_fem(1:nln,:);
end

if n_subdom == 1
    subdom{1} = [1:size(vertices,2)]';
    subdom_noOverlap{1} = 1:size(vertices,2);
    return;
end

%% compute adjacency matrix
A   = compute_adjacency(vertices, elements, dim);

%% geometric partitioning without overlap
map = metis_to_matlab(A,n_subdom,1);

%% add overlap to the partition

% subdomains numbered from 1 to n_subdom
map = map + 1; 

parfor i = 1 : n_subdom
    subdom_noOverlap{i} = find(map == i);
end

% find to which elements each node belongs
noe = size(elements,2);

nln = 3 + (dim-2);

[A_elem, ~, ~, A_elemC] = compute_adjacency_elements(vertices,elements,dim);


elem_iD          = zeros(1,noe);
elem_iD_aux      = zeros(1,noe);
elem_overlap_iD  = zeros(1,noe);
elem_interior_iD = zeros(1,noe);

subdom_elem_overlap = cell(1,n_subdom);
for ie = 1 : noe
    dof                     =  elements(1:nln,ie);
    tmp_overlap = 1;
    for kk = 2 : nln
        if map(dof(kk)) ~= map(dof(kk-1))
            tmp_overlap = 0;
        end
    end
    if tmp_overlap == 0
        for kk = 1 : nln
            subdom_elem_overlap{map(dof(kk))} = [subdom_elem_overlap{map(dof(kk))} ie];
        end
    end
    elem_iD_aux(1,ie)       =  tmp_overlap;
    elem_iD(1,ie)           =  sum(map(dof));
end
indx_overlap  = find(elem_iD_aux==0);
indx_interior = find(elem_iD_aux==1);

elem_overlap_iD(indx_overlap)   = 1;
elem_interior_iD(indx_interior) = elem_iD(indx_interior);
elem_interior_iD                = elem_interior_iD / nln;


if visual
    
    warning off
    mkdir(folder)
    warning on
    exportSubDomains(dim, elem_iD, elem_overlap_iD, elem_interior_iD, vertices,elements,[folder,'/Domain_Decomposition_overlap1'])

end



parfor i = 1 : n_subdom
    subdom_elem_interior{i} = find(elem_interior_iD == i);
end

subdom = cell(1,n_subdom);
for i = 1 : n_subdom
    dof_i     = elements_fem(:,subdom_elem_interior{i});
    dof_o     = elements_fem(:,subdom_elem_overlap{i});
    
    subdom{i} = unique([dof_i(:); dof_o(:)]);
end


if overlap > 1
    
    % loop over overlap level
    for j = 2 : overlap
        
        for i = 1 : n_subdom
            % find neighbors to a given subdomain
            [row, col]  = find(A_elem(subdom_elem_overlap{i},:));
            list_neighb = col;
            list_neighb = unique(list_neighb);
            
            % add neighbours
            subdom_elem_overlap{i}   = [subdom_elem_overlap{i} list_neighb'];
            subdom_elem_overlap{i}   = unique(subdom_elem_overlap{i});
            
            % update nodal subdomains
            dof_i     = elements_fem(:,subdom_elem_interior{i});
            dof_o     = elements_fem(:,subdom_elem_overlap{i});
            
            subdom{i} = unique([dof_i(:); dof_o(:)]);
                        
        end
        if visual
            elem_overlap_iD  = zeros(1,noe);
            for  i = 1 : n_subdom
                indx_overlap  = subdom_elem_overlap{i};
                elem_overlap_iD(indx_overlap)   = 1;
            end
            indx_overlap = find(elem_overlap_iD==1);
            elem_interior_iD(indx_overlap) = 0;
            exportSubDomains(dim, elem_iD, elem_overlap_iD, elem_interior_iD, vertices,elements,[folder,'/Domain_Decomposition_overlap',num2str(j)])
            
        end
        
    end
end
