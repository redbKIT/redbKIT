function [IDEIM_elem, DEIM_nodes, IDEIM_boundary] = CSM_DEIM_Index_to_Elements(flag, IDEIM, ndf, node_to_element, node_to_boundary, internal_vertices, nov, dim)
%CSM_DEIM_INDEX_TO_ELEMENTS detect the elements to which the nodes seletected 
%by the DEIM algorithm belong 

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

dim = dim; % to take into account for the pressure dofs as well

switch flag
      
      case 'rhs'
            
            tmp               = sparse(ndf,1);
            tmp(IDEIM)        = 1;
            
            tmp2              = sparse(dim*nov,1);
            tmp2(internal_vertices) = tmp;
            
            indx = [1:nov];
            row  = [];
            
            for k = 1 : dim
                  row  = [row; find(tmp2(indx))];
                  indx = indx + nov;
            end
            
            DEIM_nodes        = unique([row]);
            
      case 'rhs2'
            
            tmp2              = sparse(dim*nov,1);
            tmp2(IDEIM)       = 1;
            
            indx = [1:nov];
            row  = [];
            
            for k = 1 : dim
                  row  = [row; find(tmp2(indx))];
                  indx = indx + nov;
            end
            
            DEIM_nodes        = unique([row]);
            
      case 'matrix'
            
            tmp               = sparse(ndf*ndf,1);
            tmp(IDEIM)        = 1;
            
            
            tmp2              = sparse(dim*nov,dim*nov);
            tmp2(internal_vertices,internal_vertices) = reshape(tmp,ndf,ndf);
            
            row   = [];
            col   = [];
            indxR = [1:nov];            
            
            for k1 = 1 : dim
                  
                  indxC = [1:nov];
                  
                  for k2 = 1 : dim
                          
                         [row_tmp, col_tmp]    = find(tmp2(indxR,indxC));
                         
                         row     = [row; row_tmp ];
                         col     = [col; col_tmp ];  
                         
                         indxC = indxC + nov;
                  end
                  indxR = indxR + nov;
            end
            
            DEIM_nodes        = unique([row col]); %% PLEASE CHECK concatenation dimension
            
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

