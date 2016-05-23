function [A_in, F_in, u_Dirichlet] =  ADR_ApplyBC(A, F, FE_SPACE, MESH, DATA, t)
%ADR_APPLYBC apply boundary conditions for Advc-Diff-React problem in 2D/3D
%
%   [A_IN, F_IN, U_DIRICHLET] = ADR_APPLYBC(A, F, MESH, DATA) given an
%   assembled matrix A, righ-hand side vector F, a MESH data structure and
%   a DATA structure, applies Neumann, Robin and Dirichlet boundary
%   conditions. It returns the matrix A_IN (matrix A + BCs then restricted
%   to internal dofs), the vector F_IN (vector F + BCs then restricted
%   to internal dofs) and the vector U_DIRICHLET containing the
%   Dirichlet datum evaluated in the Dirichlet dofs.
%
%   Warning: in case of assembling two vectors corresponding to Neumann and
%   Dirichlet data, respectively, it is preferable to process first Neumann
%   data and then set Data.bcNeu to zero before processing Dirichlet data.
%
%   [A_IN, F_IN, U_DIRICHLET] = ADR_ApplyBC(A, F, MESH, DATA, T) as
%   before, but with the additional input T (time) for time-dependent
%   problems.

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

if nargin < 6
    t = [];
end

if isempty(A)
    A = sparse(MESH.numNodes, MESH.numNodes);
end

if isempty(F)
    F = sparse(MESH.numNodes, 1);
end

param = DATA.param;

switch MESH.dim
    case 2
        %% Neumann condition
        if ~isempty(MESH.Neumann_side)
            
            [csi,wi]       =  xwgl(FE_SPACE.quad_order, 0, 1);
            [phi]          =  fem_basis(MESH.dim, FE_SPACE.fem, [csi; 0*csi], 1);
            eta            =  1 - csi;
            nqn            =  length(csi);
            
            nof         = length(MESH.Neumann_side);
            nbn         = MESH.numBoundaryDof;
            
            Rrows       = zeros(nbn*nof,1);
            Rcoef       = Rrows;
            
            xlt = zeros(nof,nqn); ylt = xlt;
            coord_ref = [eta; csi];
            for j = 1 : 2
                dof = MESH.boundaries(j,MESH.Neumann_side);
                vtemp = MESH.vertices(1,dof);
                xlt = xlt + vtemp'*coord_ref(j,:);
                vtemp = MESH.vertices(2,dof);
                ylt = ylt + vtemp'*coord_ref(j,:);
            end
            
            u_Neumann = DATA.bcNeu(xlt,ylt,t,param);
            one       = ones(nof,nqn);
            u_Neumann = u_Neumann.*one;
            
            x    =  MESH.vertices(1,MESH.boundaries(1:2, MESH.Neumann_side));
            y    =  MESH.vertices(2,MESH.boundaries(1:2, MESH.Neumann_side));
            
            side_length = sqrt((x(2:2:end)-x(1:2:end-1)).^2+(y(2:2:end)-y(1:2:end-1)).^2);
            
            for l = 1 : nof
                face = MESH.Neumann_side(l);
                
                u_Neumann_loc  = u_Neumann(l,:).*wi;
                u_Neumann_loc  = u_Neumann_loc(1,:)';
                
                Rrows(1+(l-1)*nbn:l*nbn)    = MESH.boundaries(1:nbn,face);
                Rcoef(1+(l-1)*nbn:l*nbn)    = side_length(l)*phi*u_Neumann_loc;
            end
            F = F + sparse(Rrows,1,Rcoef,MESH.numNodes,1);
            
        end
        
        %% Robin condition
        if ~isempty(MESH.Robin_side)
            
            [csi,wi]       =  xwgl(FE_SPACE.quad_order, 0, 1);
            [phi]          =  fem_basis(MESH.dim, FE_SPACE.fem, [csi;0*csi], 1);
            eta            =  1 - csi;
            nqn            =  length(csi);
            
            
            nof         = length(MESH.Robin_side);
            nbn         = MESH.numBoundaryDof;
            
            Arows       = zeros(nbn*nbn*nof,1);
            Acols       = Arows;
            Acoef       = Arows;
            Rrows       = zeros(nbn*nof,1);
            Rcoef       = Rrows;
            [rows,cols] = meshgrid(1:nbn,1:nbn);
            rows        = rows(:);
            cols        = cols(:);
            
            xlt = zeros(nof,nqn); ylt = xlt;
            coord_ref = [eta; csi];
            for j = 1 : 2
                dof = MESH.boundaries(j,MESH.Robin_side);
                vtemp = MESH.vertices(1,dof);
                xlt = xlt + vtemp'*coord_ref(j,:);
                vtemp = MESH.vertices(2,dof);
                ylt = ylt + vtemp'*coord_ref(j,:);
            end
            
            u_Robin = DATA.bcRob_fun(xlt,ylt,t,param);
            alphaR  = DATA.bcRob_alpha(xlt,ylt,t,param);
            
            one       = ones(nof,nqn);
            u_Robin   = u_Robin.*one;
            alphaR    = alphaR.*one;
            
            x    =  MESH.vertices(1,MESH.boundaries(1:2, MESH.Robin_side));
            y    =  MESH.vertices(2,MESH.boundaries(1:2, MESH.Robin_side));
            
            side_length = sqrt((x(2:2:end)-x(1:2:end-1)).^2+(y(2:2:end)-y(1:2:end-1)).^2);
            
            for l = 1 : nof
                face = MESH.Robin_side(l);
                
                u_Robin_loc  = u_Robin(l,:).*wi;
                u_Robin_loc  = u_Robin_loc(1,:).';
                
                alphaR_loc   = alphaR(l,:).*wi;
                alphaR_loc   = alphaR_loc(ones(nbn, 1),:);
                
                Rrows(1+(l-1)*nbn:l*nbn)    = MESH.boundaries(1:nbn,face);
                Rcoef(1+(l-1)*nbn:l*nbn)    = side_length(l)*phi*u_Robin_loc;
                
                Arows(1+(l-1)*nbn*nbn:l*nbn*nbn)   =  MESH.boundaries(rows,face);
                Acols(1+(l-1)*nbn*nbn:l*nbn*nbn)   =  MESH.boundaries(cols,face);
                Acoef(1+(l-1)*nbn*nbn:l*nbn*nbn)   =  side_length(l)*(alphaR_loc.*phi)*phi';
                
            end
            
            A = A + sparse(Arows,Acols,Acoef,MESH.numNodes,MESH.numNodes);
            F = F + sparse(Rrows,1,Rcoef,MESH.numNodes,1);
            
        end
        
        
        %% Dirichlet condition
        if ~isempty(MESH.Dirichlet_dof)
            
            x           = MESH.nodes(1,MESH.Dirichlet_dof);
            y           = MESH.nodes(2,MESH.Dirichlet_dof);
            u_Dirichlet = DATA.bcDir(x,y,t,param);
            
            F_in = F(MESH.internal_dof)...
                -A(MESH.internal_dof,MESH.Dirichlet_dof)*u_Dirichlet.';
            
            A_in = A(MESH.internal_dof,MESH.internal_dof);
            
        else
            u_Dirichlet        = [];
            F_in               = F;
            A_in               = A;
        end
        
    case 3
        
        %% Neumann condition
        if ~isempty(MESH.Neumann_side)
            [quad_points, wi] = quadrature(MESH.dim-1, FE_SPACE.quad_order);
            csi = quad_points(1,:);
            eta = quad_points(2,:);
            [phi]          =  fem_basis(MESH.dim, FE_SPACE.fem, [csi; eta; 0*eta], 1);
            eta1           =  1-csi-eta;
            nqn            =  length(wi);
            
            nof         = length(MESH.Neumann_side);
            nbn         = MESH.numBoundaryDof;
            
            Rrows       = zeros(nbn*nof,1);
            Rcoef       = Rrows;
            
            xlt = zeros(nof,nqn); ylt = xlt; zlt = xlt;
            coord_ref = [eta1; csi; eta];
            for j = 1 : 2
                dof = MESH.boundaries(j,MESH.Neumann_side);
                vtemp = MESH.vertices(1,dof);
                xlt = xlt + vtemp'*coord_ref(j,:);
                vtemp = MESH.vertices(2,dof);
                ylt = ylt + vtemp'*coord_ref(j,:);
                vtemp = MESH.vertices(3,dof);
                zlt = zlt + vtemp'*coord_ref(j,:);
            end
            
            u_Neumann = DATA.bcNeu(xlt,ylt,zlt,t,param);
            one       = ones(nof,nqn);
            u_Neumann = u_Neumann.*one;
            
            x    =  MESH.vertices(1,MESH.boundaries(1:3, MESH.Neumann_side));
            y    =  MESH.vertices(2,MESH.boundaries(1:3, MESH.Neumann_side));
            z    =  MESH.vertices(3,MESH.boundaries(1:3, MESH.Neumann_side));
            
            areav = cross(  [x(2:3:end)-x(1:3:end);  y(2:3:end)-y(1:3:end);  z(2:3:end)-z(1:3:end)], ...
                [x(3:3:end)-x(1:3:end);  y(3:3:end)-y(1:3:end);  z(3:3:end)-z(1:3:end)]);
            
            for l = 1 : nof
                
                area   = 0.5*norm(areav(:,l));
                detjac = 2*area;
                
                face = MESH.Neumann_side(l);
                
                u_Neumann_loc  = u_Neumann(l,:).*wi;
                u_Neumann_loc  = u_Neumann_loc(1,:)';
                
                Rrows(1+(l-1)*nbn:l*nbn)    = MESH.boundaries(1:nbn,face);
                Rcoef(1+(l-1)*nbn:l*nbn)    = detjac*phi*u_Neumann_loc;
            end
            F = F + sparse(Rrows,1,Rcoef,MESH.numNodes,1);
        end
        
        %% Robin condition
        if ~isempty(MESH.Robin_side)
            error('3D robin BC to be implemented')
        end
        
        
        %% Dirichlet condition
        if ~isempty(MESH.Dirichlet_dof)
            
            x           = MESH.nodes(1,MESH.Dirichlet_dof);
            y           = MESH.nodes(2,MESH.Dirichlet_dof);
            z           = MESH.nodes(3,MESH.Dirichlet_dof);
            u_Dirichlet = DATA.bcDir(x,y,z,t,param);
            
            F_in = F(MESH.internal_dof)...
                -A(MESH.internal_dof,MESH.Dirichlet_dof)*u_Dirichlet.';
            
            A_in = A(MESH.internal_dof,MESH.internal_dof);
            
        else
            u_Dirichlet        = [];
            F_in               = F;
            A_in               = A;
        end
        
end


end
