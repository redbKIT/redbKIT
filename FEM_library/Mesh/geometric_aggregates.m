function [subdom_noOverlap] = geometric_aggregates(A_elem, vertices, elements, dim, n_subdom, elements_fem, out_filename)
%GEOMETRIC_DOMAIN_DECOMPOSITION builds mesh subdomains with and without overlap 
%using Metis Library
%
%   [subdom, subdom_noOverlap] = GEOMETRIC_DOMAIN_DECOMPOSITION(vertices, elements, dim, n_subdom, overlap, visual, folder, elements_fem)
%
%   see also metis_to_matlab.m

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

if n_subdom == 1
    subdom_noOverlap{1} = 1:size(nodes,2);
    return;
end

if nargin < 7
    out_filename = 'Aggregates';
end

% if nargin < 7
%     fem = 'P1';
% end
% 
% %% compute adjacency matrix
% if isempty( A )
%     A   = compute_adjacency(nodes, elements, dim, fem);
% end
% 
% %% geometric partitioning without overlap
% map = metis_to_matlab(A, n_subdom, 1);
% 
% %% add overlap to the partition
% 
% % subdomains numbered from 1 to n_subdom
% map = map + 1; 
% 
% parfor i = 1 : n_subdom
%     subdom_noOverlap{i} = find(map == i);
% end


if isempty( A_elem )
    [~,~,~,A_elem] = compute_adjacency_elements(vertices,elements,dim);
end

mapElem = metis_to_matlab(A_elem,n_subdom,1);
mapElem = mapElem + 1;

switch dim
    case 2
        exportData=struct('iteration', {-1},...
            'vertices', {vertices'},...
            'elements', {elements(1:3,:)'},...
            'outputFile', {out_filename},...
            'title', {'Aggregates'},...
            'variableName',{{'Aggregates'}},...
            'variableType',{{'SCALARS'}},...
            'variableData',{{mapElem}});
        
        exporter2dVTK_cell(exportData);
    case 3
        exportData=struct('iteration', {-1},...
            'vertices', {vertices'},...
            'elements', {elements(1:4,:)'},...
            'outputFile', {out_filename},...
            'title', {'Aggregates'},...
            'variableName',{{'Aggregates'}},...
            'variableType',{{'SCALARS'}},...
            'variableData',{{mapElem}});
        
        exporter3dVTK_cell(exportData);
end


for i = 1 : n_subdom
    sub_elem_index = find(mapElem == i);
    sub_vertices = elements_fem(:, sub_elem_index);
    subdom_noOverlap{i} = unique( sub_vertices(:) );
    if i >  1
        for j = 1 : i-1
            subdom_noOverlap{i} = setdiff( subdom_noOverlap{i}, subdom_noOverlap{j} );
        end
    end
end

end
