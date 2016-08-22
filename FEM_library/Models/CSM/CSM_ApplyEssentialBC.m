function [A_in, F_in, u_D] =  CSM_ApplyEssentialBC(A, F, MESH, DATA, t, zero_Dirichlet)
%CSM_APPLYESSENTIALBC apply essential (Dirichlet) boundary conditions for CSM problem in 2D/3D
%
%   [A_IN, F_IN, U_DIRICHLET] = CSM_APPLYESSENTIALBC(A, F, MESH, DATA) given an
%   assembled matrix A, righ-hand side vector F, a FE_SPACE, a MESH data structure and
%   a DATA structure, applies essential boundary conditions. It returns the matrix A_IN
%   (matrix A restricted to internal dofs), the vector F_IN (vector F restricted
%   to internal dofs) and the vector U_DIRICHLET containing the
%   Dirichlet datum evaluated in the Dirichlet dofs.
%
%   [A_IN, F_IN, U_DIRICHLET] = CSM_APPLYBC(A, F, MESH, DATA, T) as
%   before, but with the additional input T (time) for time-dependent
%   problems.
%
%   [A_IN, F_IN, U_DIRICHLET] = CSM_APPLYBC(A, F, MESH, DATA, T, ZERO_DIRICHLET)
%   If ZERO_DIRICHLET = 1, applies homogeneous Dirichlet boundary
%   conditions (useful for Newton iterations). ZERO_DIRICHLET = 0 by
%   default.

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

if nargin < 5
    t = [];
end

if isempty(A)
    A = sparse(MESH.numNodes*MESH.dim, MESH.numNodes*MESH.dim);
end

if isempty(F)
    F = sparse(MESH.numNodes*MESH.dim, 1);
end

if nargin < 6
    zero_Dirichlet = 0;
end

param = DATA.param;

u_D = [];

switch MESH.dim
    
    case 2
        
        %% Dirichlet condition
        for k = 1 : 2
            if ~isempty(MESH.Dirichlet_dof_c{k})
                x           = MESH.nodes(1,MESH.Dirichlet_dof_c{k});
                y           = MESH.nodes(2,MESH.Dirichlet_dof_c{k});
                u_Dirichlet{k} = DATA.bcDir{k}(x,y,t,param);
                
            else
                u_Dirichlet{k}        = [];
            end
            u_D = [u_D; u_Dirichlet{k}'];
        end
        
    case 3
        
        %% Dirichlet condition
        for k = 1 : 3
            if ~isempty(MESH.Dirichlet_dof_c{k})
                
                x           = MESH.nodes(1,MESH.Dirichlet_dof_c{k});
                y           = MESH.nodes(2,MESH.Dirichlet_dof_c{k});
                z           = MESH.nodes(3,MESH.Dirichlet_dof_c{k});
                u_Dirichlet{k} = DATA.bcDir{k}(x,y,z,t,param);
                
            else
                u_Dirichlet{k}        = [];
            end
            u_D = [u_D; u_Dirichlet{k}'];
        end
        
end

u_D  = u_D * (1 - zero_Dirichlet);

if ~isempty( MESH.Dirichlet_dof )
    
    F_in = F(MESH.internal_dof) - A(MESH.internal_dof,MESH.Dirichlet_dof)*u_D;
    
    A_in = A(MESH.internal_dof,MESH.internal_dof);
    
else
    
    F_in = F(MESH.internal_dof);
    
    A_in = A(MESH.internal_dof,MESH.internal_dof);
    
end

end
