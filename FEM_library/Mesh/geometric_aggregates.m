function [subdom_noOverlap] = geometric_aggregates(A, nodes, elements, dim, n_subdom)
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

%% compute adjacency matrix
if isempty( A )
    A   = compute_adjacency(nodes, elements, dim, fem);
end

%% geometric partitioning without overlap
map = metis_to_matlab(A, n_subdom, 1);

%% add overlap to the partition

% subdomains numbered from 1 to n_subdom
map = map + 1; 

parfor i = 1 : n_subdom
    subdom_noOverlap{i} = find(map == i);
end

end
