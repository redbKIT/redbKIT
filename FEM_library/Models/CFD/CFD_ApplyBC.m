function [A_in, F_in, u_D] =  CFD_ApplyBC(A, F, FE_SPACE, FE_SPACE_p, MESH, DATA, t, zero_Dirichlet, U_k)
%CSM_APPLYBC apply boundary conditions for CSM problem in 2D/3D
%
%   [A_IN, F_IN, U_DIRICHLET] = CFD_ApplyBC(A, F, FE_SPACE, MESH, DATA) given an
%   assembled matrix A, righ-hand side vector F, a FE_SPACE, a MESH data structure and
%   a DATA structure, applies Neumann, Normal Pressure and Dirichlet boundary
%   conditions. It returns the matrix A_IN (matrix A + BCs then restricted
%   to internal dofs), the vector F_IN (vector F + BCs then restricted
%   to internal dofs) and the vector U_DIRICHLET containing the
%   Dirichlet datum evaluated in the Dirichlet dofs.
%
%   [A_IN, F_IN, U_DIRICHLET] = CFD_ApplyBC(A, F, FE_SPACE, MESH, DATA, T) as
%   before, but with the additional input T (time) for time-dependent
%   problems.
%
%   [A_IN, F_IN, U_DIRICHLET] = CFD_ApplyBC(A, F, FE_SPACE, MESH, DATA, T, ZERO_DIRICHLET)
%   If ZERO_DIRICHLET = 1, applies homogeneous Dirichlet boundary
%   conditions (useful for Newton iterations). ZERO_DIRICHLET = 0 by
%   default.

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

if nargin < 7
    t = [];
end

if isempty(A)
    A = sparse(FE_SPACE.numDof + FE_SPACE_p.numDof, FE_SPACE.numDof + FE_SPACE_p.numDof);
end

if isempty(F)
    F = sparse(FE_SPACE.numDof + FE_SPACE_p.numDof, 1);
end

if nargin < 8
    zero_Dirichlet = 0;
end


if nargin < 9
    U_k = zeros(FE_SPACE.numDof + FE_SPACE_p.numDof,1);
end

param = DATA.param;

u_D = [];

switch MESH.dim
    case 2
        %% Pressure condition
        for k = 1 : MESH.dim
            if ~isempty(MESH.Pressure_side{k})
                
                [csi,wi]       =  xwgl(FE_SPACE.quad_order, 0, 1);
                [phi]          =  fem_basis(MESH.dim, FE_SPACE.fem, [csi; 0*csi], 1);
                eta            =  1 - csi;
                nqn            =  length(csi);
                
                nof         = length(MESH.Pressure_side{k});
                nbn         = MESH.numBoundaryDof;
                
                Rrows       = zeros(nbn*nof,1);
                Rcoef       = Rrows;
                
                xlt = zeros(nof,nqn); ylt = xlt;
                coord_ref = [eta; csi];
                for j = 1 : 2
                    dof = MESH.boundaries(j,MESH.Pressure_side{k});
                    vtemp = MESH.vertices(1,dof);
                    xlt = xlt + vtemp'*coord_ref(j,:);
                    vtemp = MESH.vertices(2,dof);
                    ylt = ylt + vtemp'*coord_ref(j,:);
                end
                
                pressure = DATA.bcPrex(xlt,ylt,t,param);
                one       = ones(nof,nqn);
                pressure = pressure.*one;
                
                x    =  MESH.vertices(1,MESH.boundaries(1:2, MESH.Pressure_side{k}));
                y    =  MESH.vertices(2,MESH.boundaries(1:2, MESH.Pressure_side{k}));
                
                side_length = sqrt((x(2:2:end)-x(1:2:end-1)).^2+(y(2:2:end)-y(1:2:end-1)).^2);
                
                for l = 1 : nof
                    face = MESH.Pressure_side{k}(l);
                    
                    pressure_loc  = pressure(l,:).*wi;
                    pressure_loc  = pressure_loc(1,:)';
                    
                    Rrows(1+(l-1)*nbn:l*nbn)    = MESH.boundaries(1:nbn,face);
                    Rcoef(1+(l-1)*nbn:l*nbn)    = MESH.Normal_Faces(k,face)*side_length(l)*phi*pressure_loc;
                end
                F = F + sparse(Rrows+(k-1)*MESH.numNodes,1,Rcoef,FE_SPACE.numDof + FE_SPACE_p.numDof,1);
                
            end
        end
        
        %% Dirichlet condition
        if size(DATA.bcDir,1) == 1
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
        
        else
            v_D    = zeros(FE_SPACE.numDof,1);
            % loop over components
            for k = 1 : 2
                if ~isempty(MESH.Dirichlet_dof_c{k})
                    
                    % loop over dirichlet boundaries
                    for f = 1 : length(DATA.flag_dirichlet{k})
                        
                        x           = MESH.nodes(1,MESH.DiriDof_CompFlag{k,f});
                        y           = MESH.nodes(2,MESH.DiriDof_CompFlag{k,f});
                        
                        u_Dirichlet = DATA.bcDir{k,DATA.flag_dirichlet{k}(f)}(x,y,t,param);
                        
                        v_D(MESH.DiriDof_CompFlag{k,f}+(k-1)*FE_SPACE.numDofScalar) = u_Dirichlet';
                        
                    end
                end
            end
            u_D = v_D(MESH.Dirichlet_dof);
        end
        
        
    case 3
        
        %% Neumann condition
        for k = 1 : MESH.dim
            if ~isempty(MESH.Neumann_side{k})
                
                [quad_points, wi] = quadrature(MESH.dim-1, FE_SPACE.quad_order);
                csi = quad_points(1,:);
                eta = quad_points(2,:);
                [phi]          =  fem_basis(MESH.dim, FE_SPACE.fem, [csi; eta; 0*eta], 1);
                eta1           =  1-csi-eta;
                nqn            =  length(wi);
                
                nof         = length(MESH.Neumann_side{k});
                nbn         = MESH.numBoundaryDof;
                
                Rrows       = zeros(nbn*nof,1);
                Rcoef       = Rrows;
                
                xlt = zeros(nof,nqn); ylt = xlt; zlt = xlt;
                coord_ref = [eta1; csi; eta];
                for j = 1 : 2
                    dof = MESH.boundaries(j,MESH.Neumann_side{k});
                    vtemp = MESH.vertices(1,dof);
                    xlt = xlt + vtemp'*coord_ref(j,:);
                    vtemp = MESH.vertices(2,dof);
                    ylt = ylt + vtemp'*coord_ref(j,:);
                    vtemp = MESH.vertices(3,dof);
                    zlt = zlt + vtemp'*coord_ref(j,:);
                end
                
                u_Neumann = DATA.bcNeu{k}(xlt,ylt,zlt,t,param);
                one       = ones(nof,nqn);
                u_Neumann = u_Neumann.*one;
                
                x    =  MESH.vertices(1,MESH.boundaries(1:3, MESH.Neumann_side{k}));
                y    =  MESH.vertices(2,MESH.boundaries(1:3, MESH.Neumann_side{k}));
                z    =  MESH.vertices(3,MESH.boundaries(1:3, MESH.Neumann_side{k}));
                
                areav = cross(  [x(2:3:end)-x(1:3:end);  y(2:3:end)-y(1:3:end);  z(2:3:end)-z(1:3:end)], ...
                  [x(3:3:end)-x(1:3:end);  y(3:3:end)-y(1:3:end);  z(3:3:end)-z(1:3:end)]);
                
                for l = 1 : nof
                    
                    area   = 0.5*norm(areav(:,l));
                    detjac = 2*area;
                    
                    face = MESH.Neumann_side{k}(l);
                    
                    u_Neumann_loc  = u_Neumann(l,:).*wi;
                    u_Neumann_loc  = u_Neumann_loc(1,:)';
                    
                    Rrows(1+(l-1)*nbn:l*nbn)    = MESH.boundaries(1:nbn,face);
                    Rcoef(1+(l-1)*nbn:l*nbn)    = detjac*phi*u_Neumann_loc;
                end
                F = F + sparse(Rrows+(k-1)*MESH.numNodes,1,Rcoef,FE_SPACE.numDof + FE_SPACE_p.numDof,1);
                
            end
        end
        
        
        %% Pressure condition
        for k = 1 : MESH.dim
            if ~isempty(MESH.Pressure_side{k})
                
                [quad_points, wi] = quadrature(MESH.dim-1, FE_SPACE.quad_order);
                csi = quad_points(1,:);
                eta = quad_points(2,:);
                [phi]          =  fem_basis(MESH.dim, FE_SPACE.fem, [csi; eta; 0*eta], 1);
                eta1           =  1-csi-eta;
                nqn            =  length(wi);
                
                nbn         = MESH.numBoundaryDof;
                
                for flag = 1 : length(DATA.flag_pressure{k})
                    
                    nof         = length(MESH.Pressure_side_CompFlag{k,flag});

                    Rrows       = zeros(nbn*nof,1);
                    Rcoef       = Rrows;
                    
                    xlt = zeros(nof,nqn); ylt = xlt; zlt = xlt;
                    coord_ref = [eta1; csi; eta];
                    for j = 1 : 2
                        dof = MESH.boundaries(j,MESH.Pressure_side_CompFlag{k,flag});
                        vtemp = MESH.vertices(1,dof);
                        xlt = xlt + vtemp'*coord_ref(j,:);
                        vtemp = MESH.vertices(2,dof);
                        ylt = ylt + vtemp'*coord_ref(j,:);
                        vtemp = MESH.vertices(3,dof);
                        zlt = zlt + vtemp'*coord_ref(j,:);
                    end
                    
                    if length(DATA.bcPrex) == 1
                        pressure = DATA.bcPrex(xlt,ylt,zlt,t,param);
                    else
                        pressure = DATA.bcPrex{DATA.flag_pressure{k}(flag)}(xlt,ylt,zlt,t,param);
                    end
                    
                    one      = ones(nof,nqn);
                    pressure = pressure.*one;
                    
                    x    =  MESH.vertices(1,MESH.boundaries(1:3, MESH.Pressure_side_CompFlag{k,flag}));
                    y    =  MESH.vertices(2,MESH.boundaries(1:3, MESH.Pressure_side_CompFlag{k,flag}));
                    z    =  MESH.vertices(3,MESH.boundaries(1:3, MESH.Pressure_side_CompFlag{k,flag}));
                    
                    areav = cross(  [x(2:3:end)-x(1:3:end);  y(2:3:end)-y(1:3:end);  z(2:3:end)-z(1:3:end)], ...
                        [x(3:3:end)-x(1:3:end);  y(3:3:end)-y(1:3:end);  z(3:3:end)-z(1:3:end)]);
                    
                    for l = 1 : nof
                        
                        area   = 0.5*norm(areav(:,l));
                        detjac = 2*area;
                        
                        face = MESH.Pressure_side_CompFlag{k,flag}(l);
                        
                        pressure_loc  = pressure(l,:).*wi;
                        pressure_loc  = pressure_loc(1,:)';
                        
                        Rrows(1+(l-1)*nbn:l*nbn)    = MESH.boundaries(1:nbn,face);
                        Rcoef(1+(l-1)*nbn:l*nbn)    = MESH.Normal_Faces(k,face)*detjac*phi*pressure_loc;
                    end
                    F = F + GlobalAssemble(Rrows+(k-1)*MESH.numNodes,1,Rcoef,FE_SPACE.numDof + FE_SPACE_p.numDof,1);
                    
                end
            end
        end
        
        %% Resistance condition
        if ~isempty(MESH.Resistance_side)
            
            [quad_points, wi] = quadrature(MESH.dim-1, FE_SPACE.quad_order);
            csi = quad_points(1,:);
            eta = quad_points(2,:);
            [phi]    =  fem_basis(MESH.dim, FE_SPACE.fem, [csi; eta; 0*eta], 1);
            nbn      = MESH.numBoundaryDof;
                        
            %fprintf('\n');
            for flag = 1 : length(DATA.flag_resistance)
                
                [FlowRate, Area, P_average]    = CFD_computeFlowRate(U_k, MESH, FE_SPACE, FE_SPACE_p, DATA.flag_resistance(flag));
                nof         = length(MESH.Resistance_side_Flag{flag});
                %fprintf('\n Flag %d, Area  = %2.6f', DATA.flag_resistance(flag), Area);
                
                if isfield(DATA, 'OutFlow_ComputeApproximateResistance')
                    if DATA.OutFlow_ComputeApproximateResistance == true
                        Radius   = sqrt( Area / pi);
                        beta     = sqrt( pi ) * DATA.OutFlow_WallThicknessPercentage * Radius * DATA.OutFlow_Young / (1 - DATA.OutFlow_Poisson^2);
                        DATA.OutFlow_Resistance(flag) = sqrt(DATA.density * beta) / sqrt(2) / Area^(5/4);
                    end
                end
                
                Rrows       = zeros(nbn*nof,1);
                Rcoef       = Rrows;
                
                pressure = DATA.time.dt / (DATA.time.dt + DATA.OutFlow_Resistance(flag) * DATA.OutFlow_Capacity(flag)) ...
                                 * ( FlowRate * DATA.OutFlow_Resistance(flag) + P_average * DATA.OutFlow_Resistance(flag) * DATA.OutFlow_Capacity(flag) / DATA.time.dt ...
                                     + DATA.OutFlow_RefPressure(t) );
                
                x    =  MESH.vertices(1,MESH.boundaries(1:3, MESH.Resistance_side_Flag{flag}));
                y    =  MESH.vertices(2,MESH.boundaries(1:3, MESH.Resistance_side_Flag{flag}));
                z    =  MESH.vertices(3,MESH.boundaries(1:3, MESH.Resistance_side_Flag{flag}));
                
                areav = cross(  [x(2:3:end)-x(1:3:end);  y(2:3:end)-y(1:3:end);  z(2:3:end)-z(1:3:end)], ...
                    [x(3:3:end)-x(1:3:end);  y(3:3:end)-y(1:3:end);  z(3:3:end)-z(1:3:end)]);
                
                for k = 1 : MESH.dim
                    
                    for l = 1 : nof
                        
                        area   = 0.5*norm(areav(:,l));
                        detjac = 2*area;
                        
                        face = MESH.Resistance_side_Flag{flag}(l);
                        
                        pressure_loc  = pressure * wi';
                        
                        Rrows(1+(l-1)*nbn:l*nbn)    = MESH.boundaries(1:nbn,face);
                        Rcoef(1+(l-1)*nbn:l*nbn)    = - MESH.Normal_Faces(k,face)*detjac*phi*pressure_loc;
                    end
                    F = F + GlobalAssemble(Rrows+(k-1)*MESH.numNodes,1,Rcoef,FE_SPACE.numDof + FE_SPACE_p.numDof,1);
                end
                
            end
            
        end
        
        %% Absorbing condition
        if ~isempty(MESH.Absorbing_side)
            
            [quad_points, wi] = quadrature(MESH.dim-1, FE_SPACE.quad_order);
            csi = quad_points(1,:);
            eta = quad_points(2,:);
            [phi]    =  fem_basis(MESH.dim, FE_SPACE.fem, [csi; eta; 0*eta], 1);
            nbn      =  MESH.numBoundaryDof;
            
            %fprintf('\n');
            for flag = 1 : length(DATA.flag_absorbing)
                
                [FlowRate, Area]    = CFD_computeFlowRate(U_k, MESH, FE_SPACE, FE_SPACE_p, DATA.flag_absorbing(flag));
                nof         = length(MESH.Absorbing_side_Flag{flag});
                
                Rrows       = zeros(nbn*nof,1);
                Rcoef       = Rrows;
                
                Radius   = sqrt(Area / pi);
                beta     = sqrt( pi ) * DATA.OutFlow_WallThicknessPercentage * Radius * DATA.OutFlow_Young / (1 - DATA.OutFlow_Poisson^2);
                pressure = sqrt(DATA.density * beta) / sqrt(2) / Area^(5/4) * FlowRate + DATA.OutFlow_RefPressure(t);
                % pressure = ( sqrt(DATA.density) / (2*sqrt(2)) * FlowRate / Area + sqrt( beta * sqrt(Area) ) )^2  - beta * sqrt(Area) + DATA.OutFlow_RefPressure(t);
                %fprintf('\n Flag %d, beta = %2.6f, Equivalent Resistance  = %2.6f',  DATA.flag_absorbing(flag), beta,  sqrt(DATA.density * beta) / sqrt(2) / Area^(5/4) ) ;

                x    =  MESH.vertices(1,MESH.boundaries(1:3, MESH.Absorbing_side_Flag{flag}));
                y    =  MESH.vertices(2,MESH.boundaries(1:3, MESH.Absorbing_side_Flag{flag}));
                z    =  MESH.vertices(3,MESH.boundaries(1:3, MESH.Absorbing_side_Flag{flag}));
                
                areav = cross(  [x(2:3:end)-x(1:3:end);  y(2:3:end)-y(1:3:end);  z(2:3:end)-z(1:3:end)], ...
                    [x(3:3:end)-x(1:3:end);  y(3:3:end)-y(1:3:end);  z(3:3:end)-z(1:3:end)]);
                
                for k = 1 : MESH.dim
                    
                    for l = 1 : nof
                        
                        area   = 0.5*norm(areav(:,l));
                        detjac = 2*area;
                        
                        face = MESH.Absorbing_side_Flag{flag}(l);
                        
                        pressure_loc  = pressure * wi';
                        
                        Rrows(1+(l-1)*nbn:l*nbn)    = MESH.boundaries(1:nbn,face);
                        Rcoef(1+(l-1)*nbn:l*nbn)    = - MESH.Normal_Faces(k,face)*detjac*phi*pressure_loc;
                    end
                    F = F + GlobalAssemble(Rrows+(k-1)*MESH.numNodes,1,Rcoef,FE_SPACE.numDof + FE_SPACE_p.numDof,1);
                end
                
            end
            
        end
        
        %% Dirichlet condition
        if size(DATA.bcDir,1) == 1
            
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
            
        else
            
            v_D    = zeros(FE_SPACE.numDof,1);
            % loop over components
            for k = 1 : 3
                if ~isempty(MESH.Dirichlet_dof_c{k})
                    
                    % loop over dirichlet boundaries
                    for f = 1 : length(DATA.flag_dirichlet{k})
                        
                        x           = MESH.nodes(1,MESH.DiriDof_CompFlag{k,f});
                        y           = MESH.nodes(2,MESH.DiriDof_CompFlag{k,f});
                        z           = MESH.nodes(3,MESH.DiriDof_CompFlag{k,f});
                        
                        u_Dirichlet = DATA.bcDir{k,DATA.flag_dirichlet{k}(f)}(x,y,z,t,param);
                        
                        v_D(MESH.DiriDof_CompFlag{k,f}+(k-1)*FE_SPACE.numDofScalar) = u_Dirichlet';
                        
                    end
                end
                v_D(MESH.ringDofs{k}+(k-1)*FE_SPACE.numDofScalar) = 0; % set zero velocity on rings
            end
            %[Q, Area] = CFD_computeFlowRate([v_D;zeros(FE_SPACE_p.numDof,1)], MESH, FE_SPACE, FE_SPACE_p, 3)
            u_D = v_D(MESH.Dirichlet_dof);% [Q, Area] = CFD_computeFlowRate(v_D, MESH, FE_SPACE, 1);    
        end
end

u_D  = u_D * (1 - zero_Dirichlet);

F_in = F(MESH.internal_dof) - A(MESH.internal_dof,MESH.Dirichlet_dof)*u_D;

A_in = A(MESH.internal_dof,MESH.internal_dof);

end
