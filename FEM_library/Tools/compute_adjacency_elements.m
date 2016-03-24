function [ A, node_to_element, node_to_boundaries ] = compute_adjacency_elements(vertices, elements, dim, boundaries)
%COMPUTE_ADJACENCY compute adjacency matrix for 2d or 3d TRI/TET P1 mesh

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

noe = size(elements,2);
nov = size(vertices,2);
node_to_element = cell(1,nov);

nln = 3 + (dim-2);

for ie = 1 : noe
    dof = elements(:,ie);
    
    for k = 1 : nln
       node_to_element{dof(k)} = [node_to_element{dof(k)} ie];
    end    
end

if nargin ==4 && ~isempty(boundaries)
      
      nbp                = size(boundaries,2);
      node_to_boundaries = cell(1,nov);
      nbn                = 2 + (dim-2);
      for ib = 1 : nbp
            
            dof = boundaries(1:nbn,ib);
            
            for k = 1 : nbn
                  node_to_boundaries{dof(k)} = [node_to_boundaries{dof(k)} ib];
            end
      end

end


row  = zeros(nln*noe,1);
col  = row;
coef = row;

ii   = 0;

for i = 1 : nov
    
    iii = [ii+1 : ii+length(node_to_element{i})^2];
        
    elem_n = node_to_element{i};
    
    row_tmp  = elem_n;
    col_tmp  = elem_n;
    
    [row_tmp,col_tmp]  =  meshgrid(row_tmp, col_tmp);

    row(iii)    = row_tmp(:);
    col(iii)    = col_tmp(:);
    
    coef(iii)   = 1;
    
    ii  = iii(end);

end

A = sparse(row,col,coef,noe,noe);

[row,col,coef] = find(A);

indx = find(coef>1);
coef(indx) = 1;

A = sparse(row,col,coef,noe,noe);

clear row col coef

end
