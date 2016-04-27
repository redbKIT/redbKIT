function [VIS, R_DD ] = STR_setup_Precon( VIS, LinSolveOpt )

%   Author: F. Negri (federico.negri@epfl.ch) 2014
%   Copyright (C) Federico Negri, CMCS, EPFL
%


switch LinSolveOpt.solver
      
      case 'DD_Schwarz'
            
            fprintf('\nDD_Schwarz preconditioner: subdomain partitioning, build restriction-prolongation operators\n')
            
            
            [VIS.subdom, R_DD] = STR_overlapping_DD(VIS.verticesP1, VIS.elementsP1, ...
                  VIS.dim, LinSolveOpt.partit_level, LinSolveOpt.overlap_level, VIS.internal_vertices, VIS.I, VIS.nln, 1, VIS.fem, VIS.vertices, VIS.elements);
            
 
      otherwise
            
            R_DD             = [];
end


end
