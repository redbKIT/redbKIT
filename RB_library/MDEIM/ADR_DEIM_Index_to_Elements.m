function [IDEIM_elem, DEIM_nodes, IDEIM_boundary] = ADR_DEIM_Index_to_Elements(flag, IDEIM, ndf, node_to_element, node_to_boundary, internal_vertices, nov)
%ADR_DEIM_INDEX_TO_ELEMENTS detect the elements to which the nodes seletected 
%by the DEIM algorithm belong 
%
%   IDEIM: indices vector provided by DEIM.m
%   NDF:   number of DOFs 
%   NODE_TO_ELEMENT: cell array given by compute_adjacency_elements.m
%   INTERNAL_VERTICES: index arrays of internal_vertices (since we use lifting to treat Dirichlet BCs)
%   NOV: number of vertices of the mesh

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 


switch flag
      
      case 'rhs'
            
            tmp               = sparse(ndf,1);
            tmp(IDEIM)        = 1;
            
            tmp2              = sparse(nov,1);
            tmp2(internal_vertices) = tmp;
            [row]    = find(tmp2);
            DEIM_nodes        = unique([row]);
            
      case 'matrix'
            
            tmp               = sparse(ndf*ndf,1);
            tmp(IDEIM)        = 1;
            
            
            tmp2              = sparse(nov,nov);
            tmp2(internal_vertices,internal_vertices) = reshape(tmp,ndf,ndf);
            [row,col,coef]    = find(tmp2);
            DEIM_nodes        = unique([row col]);
            
end



IDEIM_elem = [];
for i = 1 : length(DEIM_nodes)
      IDEIM_elem  = [IDEIM_elem node_to_element{DEIM_nodes(i)}];
end

IDEIM_elem = unique(IDEIM_elem);



IDEIM_boundary = [];
for i = 1 : length(DEIM_nodes)
      IDEIM_boundary  = [IDEIM_boundary node_to_boundary{DEIM_nodes(i)}];
end

IDEIM_boundary = unique(IDEIM_boundary);

end

