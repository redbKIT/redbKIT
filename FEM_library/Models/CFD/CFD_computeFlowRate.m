function [Q, Area, P_average] = CFD_computeFlowRate(u_n, MESH, FE_SPACE_v, FE_SPACE_p, boundary_flag)

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

nflags = length(boundary_flag);
face_list = [];
for j = 1 : nflags
    tmp            =  find(MESH.boundaries(MESH.bc_flag_row,:) == boundary_flag(j));
    face_list      =  [face_list, tmp];
end
face_list = unique(face_list);

Q         = 0;
Area      = 0;
P_average = 0;


switch MESH.dim
    
    case 3
        
        [quad_points, wi] = quadrature(MESH.dim-1, 5);
        csi = quad_points(1,:);
        eta = quad_points(2,:);
        phi            =  fem_basis(MESH.dim, FE_SPACE_v.fem, [csi; eta; 0*eta], 1);
        phi_p          =  fem_basis(MESH.dim, FE_SPACE_p.fem, [csi; eta; 0*eta], 1);

        nqn = length(wi);
        
        if ~isempty(face_list)
            
            nof         = length(face_list);
            
            x    =  MESH.nodes(1,MESH.boundaries(1:3,face_list));
            y    =  MESH.nodes(2,MESH.boundaries(1:3,face_list));
            z    =  MESH.nodes(3,MESH.boundaries(1:3,face_list));
            
            areav = cross(  [x(2:3:end)-x(1:3:end);  y(2:3:end)-y(1:3:end);  z(2:3:end)-z(1:3:end)], ...
                [x(3:3:end)-x(1:3:end);  y(3:3:end)-y(1:3:end);  z(3:3:end)-z(1:3:end)]);
            
            for l = 1 : nof
                
                face   = face_list(l);
                area   = 0.5*norm(areav(:,l));
                detjac = 2*area;
                
                dofs      = MESH.boundaries(1:MESH.numBoundaryDof,face);
                dofs_p    = MESH.boundaries(1:3,face);
                
                uhq = (u_n(dofs)'*phi);
                vhq = (u_n(dofs+FE_SPACE_v.numDofScalar)'*phi);
                whq = (u_n(dofs+2*FE_SPACE_v.numDofScalar)'*phi);
                phq = (u_n(dofs_p+3*FE_SPACE_v.numDofScalar)'*phi_p);
                
                for q = 1 : nqn
                    Q    = Q + ( uhq(q)*MESH.Normal_Faces(1,face) + vhq(q)*MESH.Normal_Faces(2,face) + whq(q)*MESH.Normal_Faces(3,face)) * wi(q) * detjac;
                    Area = Area + wi(q) * detjac;
                    P_average = P_average +  phq(q) * wi(q) * detjac;
                end
                
            end
            
        end
        
    case 2
        
        warning('Average pressure computation only available in 3D');

        [csi,wi]         =  xwgl(FE_SPACE_v.quad_order, 0, 1);
        [phi]          =  fem_basis(MESH.dim, FE_SPACE_v.fem, [csi; 0*csi], 1);
                
        nqn = length(wi);
        
        if ~isempty(face_list)
            
            nof         = length(face_list);
            
            x    =  MESH.nodes(1,MESH.boundaries(1:2,face_list));
            y    =  MESH.nodes(2,MESH.boundaries(1:2,face_list));
            
            detjac = sqrt((x(2:2:end)-x(1:2:end-1)).^2+(y(2:2:end)-y(1:2:end-1)).^2);

            for l = 1 : nof
                
                face   = face_list(l);
                
                dofs    = MESH.boundaries(1:MESH.numBoundaryDof,face);
                
                uhq = (u_n(dofs)'*phi);
                vhq = (u_n(dofs+FE_SPACE_v.numDofScalar)'*phi);
                
                for q = 1 : nqn
                    Q = Q + ( uhq(q)*MESH.Normal_Faces(1,face) + vhq(q)*MESH.Normal_Faces(2,face)) * wi(q) * detjac(l);
                    Area = Area + wi(q) * detjac;
                end
                
            end
            
        end
end

P_average = P_average / Area;

end
