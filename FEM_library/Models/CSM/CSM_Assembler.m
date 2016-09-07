%CSM_ASSEMBLER assembler class for 2D/3D Computational Solid Mechanics
% CSM_ASSEMBLER methods:
%    CSM_Assembler                - constructor
%    SetMaterialParameters        - set parameters vector
%    compute_volumetric_forces    - assemble volumetric rhs contribute 
%    compute_surface_forces       - assemble surface rhs contribute 
%    compute_external_forces      - assemble all external forces
%    compute_mass                 - assemble mass matrix
%    compute_stress               - compute stress for postprocessing
%    compute_internal_forces      - assemble vector of internal forces
%    compute_jacobian             - assemble jacobian (tangent stiffness) matrix
%
% CSM_ASSEMBLER properties:
%    M_MESH             - struct containing MESH data
%    M_DATA             - struct containing DATA information
%    M_FE_SPACE         - struct containing Finite Element space data
%    M_MaterialModel    - string containing name of the material model
%    M_MaterialParam    - vector containing material parameters

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri@epfl.ch>

classdef CSM_Assembler < handle
    
    properties (GetAccess = public, SetAccess = protected)
        M_MESH;
        M_DATA;
        M_FE_SPACE;
        M_subdomain;
        M_MaterialModel;
        M_MaterialParam;
    end
   
    methods
        
        %==========================================================================
        %% Constructor
        function obj = CSM_Assembler( MESH, DATA, FE_SPACE )
            
            obj.M_MESH      = MESH;
            obj.M_DATA      = DATA;
            obj.M_FE_SPACE  = FE_SPACE;
            obj.M_MaterialModel = DATA.Material_Model;
            obj = SetMaterialParameters(obj);
            
            if obj.M_MESH.dim == 2 
                if strcmp(obj.M_MaterialModel, 'NeoHookean')
                    error('NeoHookean material law is available only for 3D simulations.')
                end
                if strcmp(obj.M_MaterialModel, 'RaghavanVorp') 
                    error('RaghavanVorp material law is available only for 3D simulations.')
                end
            end
            
        end
        
        %==========================================================================
        %% SetMaterialParameters
        function obj = SetMaterialParameters( obj )
            
            switch obj.M_MaterialModel
                case {'Linear', 'StVenantKirchhoff', 'NeoHookean'}
                    obj.M_MaterialParam = [obj.M_DATA.Young obj.M_DATA.Poisson];
                    
                case {'RaghavanVorp'}
                    obj.M_MaterialParam = [obj.M_DATA.Alpha obj.M_DATA.Beta obj.M_DATA.Bulk];
                     
                case 'SEMMT'
                    obj.M_MaterialParam = [obj.M_DATA.Young obj.M_DATA.Poisson obj.M_DATA.Stiffening_power];
            end
            
        end
        
        %==========================================================================
        %% Compute Volumetric Forces
        function F_ext = compute_volumetric_forces( obj, t )
            
            if nargin < 2 || isempty(t)
                t = [];
            end
            
            % Computations of all quadrature nodes in the elements
            coord_ref = obj.M_MESH.chi;
            switch obj.M_MESH.dim
                
                case 2
                    
                    x = zeros(obj.M_MESH.numElem,obj.M_FE_SPACE.numQuadNodes); y = x;
                    for j = 1 : 3
                        i = obj.M_MESH.elements(j,:);
                        vtemp = obj.M_MESH.vertices(1,i);
                        x = x + vtemp'*coord_ref(j,:);
                        vtemp = obj.M_MESH.vertices(2,i);
                        y = y + vtemp'*coord_ref(j,:);
                    end
                    
                    % Evaluation of external forces in the quadrature nodes
                    for k = 1 : obj.M_MESH.dim
                        f{k}  = obj.M_DATA.force{k}(x,y,t,obj.M_DATA.param);
                    end
                    
                case 3
                    
                    x = zeros(obj.M_MESH.numElem,obj.M_FE_SPACE.numQuadNodes); y = x; z = x;
                    
                    for j = 1 : 4
                        i = obj.M_MESH.elements(j,:);
                        vtemp = obj.M_MESH.vertices(1,i);
                        x = x + vtemp'*coord_ref(j,:);
                        vtemp = obj.M_MESH.vertices(2,i);
                        y = y + vtemp'*coord_ref(j,:);
                        vtemp = obj.M_MESH.vertices(3,i);
                        z = z + vtemp'*coord_ref(j,:);
                    end
                    
                    % Evaluation of external forces in the quadrature nodes
                    for k = 1 : obj.M_MESH.dim
                        f{k}  = obj.M_DATA.force{k}(x,y,z,t,obj.M_DATA.param);
                    end
                    
            end
            % C_OMP assembly, returns matrices in sparse vector format

            F_ext = [];
            for k = 1 : obj.M_MESH.dim
                
                [rowF, coefF] = CSM_assembler_ExtForces(f{k}, obj.M_MESH.elements, obj.M_FE_SPACE.numElemDof, ...
                    obj.M_FE_SPACE.quad_weights, obj.M_MESH.jac, obj.M_FE_SPACE.phi);
                
                % Build sparse matrix and vector
                F_ext    = [F_ext; GlobalAssemble(rowF, 1, coefF, obj.M_MESH.numNodes, 1)];
                
            end
            
        end
        
        %==========================================================================
        %% Compute Surface Forces
        function F = compute_surface_forces( obj, t )
            
            if nargin < 2 || isempty(t)
                t = [];
            end
            
            F = sparse(obj.M_MESH.numNodes*obj.M_MESH.dim, 1);

            switch obj.M_MESH.dim
                
                case 2
                    % Pressure condition
                    for k = 1 : obj.M_MESH.dim
                        if ~isempty(obj.M_MESH.Pressure_side{k})
                            
                            [csi,wi]       =  xwgl(obj.M_FE_SPACE.quad_order, 0, 1);
                            [phi]          =  fem_basis(obj.M_MESH.dim, obj.M_FE_SPACE.fem, [csi; 0*csi], 1);
                            eta            =  1 - csi;
                            nqn            =  length(csi);
                            
                            nof         = length(obj.M_MESH.Pressure_side{k});
                            nbn         = obj.M_MESH.numBoundaryDof;
                            
                            Rrows       = zeros(nbn*nof,1);
                            Rcoef       = Rrows;
                            
                            xlt = zeros(nof,nqn); ylt = xlt;
                            coord_ref = [eta; csi];
                            for j = 1 : 2
                                dof = obj.M_MESH.boundaries(j,obj.M_MESH.Pressure_side{k});
                                vtemp = obj.M_MESH.vertices(1,dof);
                                xlt = xlt + vtemp'*coord_ref(j,:);
                                vtemp = obj.M_MESH.vertices(2,dof);
                                ylt = ylt + vtemp'*coord_ref(j,:);
                            end
                            
                            pressure = obj.M_DATA.bcPrex(xlt,ylt,t,obj.M_DATA.param);
                            one       = ones(nof,nqn);
                            pressure = pressure.*one;
                            
                            x    =  obj.M_MESH.vertices(1,obj.M_MESH.boundaries(1:2, obj.M_MESH.Pressure_side{k}));
                            y    =  obj.M_MESH.vertices(2,obj.M_MESH.boundaries(1:2, obj.M_MESH.Pressure_side{k}));
                            
                            side_length = sqrt((x(2:2:end)-x(1:2:end-1)).^2+(y(2:2:end)-y(1:2:end-1)).^2);
                            
                            for l = 1 : nof
                                face = obj.M_MESH.Pressure_side{k}(l);
                                
                                pressure_loc  = pressure(l,:).*wi;
                                pressure_loc  = pressure_loc(1,:)';
                                
                                Rrows(1+(l-1)*nbn:l*nbn)    = obj.M_MESH.boundaries(1:nbn,face);
                                Rcoef(1+(l-1)*nbn:l*nbn)    = obj.M_MESH.Normal_Faces(k,face)*side_length(l)*phi*pressure_loc;
                            end
                            F = F + sparse(Rrows+(k-1)*obj.M_MESH.numNodes,1,Rcoef,obj.M_MESH.dim*obj.M_MESH.numNodes,1);
                            
                        end
                    end
                    
                    
                case 3
                    % Neumann condition
                    for k = 1 : obj.M_MESH.dim
                        if ~isempty(obj.M_MESH.Neumann_side{k})
                            
                            [quad_points, wi] = quadrature(obj.M_MESH.dim-1, obj.M_FE_SPACE.quad_order);
                            csi = quad_points(1,:);
                            eta = quad_points(2,:);
                            [phi]          =  fem_basis(obj.M_MESH.dim, obj.M_FE_SPACE.fem, [csi; eta; 0*eta], 1);
                            eta1           =  1-csi-eta;
                            nqn            =  length(wi);
                            
                            nof         = length(obj.M_MESH.Neumann_side{k});
                            nbn         = obj.M_MESH.numBoundaryDof;
                            
                            Rrows       = zeros(nbn*nof,1);
                            Rcoef       = Rrows;
                            
                            xlt = zeros(nof,nqn); ylt = xlt; zlt = xlt;
                            coord_ref = [eta1; csi; eta];
                            for j = 1 : 2
                                dof = obj.M_MESH.boundaries(j,obj.M_MESH.Neumann_side{k});
                                vtemp = obj.M_MESH.vertices(1,dof);
                                xlt = xlt + vtemp'*coord_ref(j,:);
                                vtemp = obj.M_MESH.vertices(2,dof);
                                ylt = ylt + vtemp'*coord_ref(j,:);
                                vtemp = obj.M_MESH.vertices(3,dof);
                                zlt = zlt + vtemp'*coord_ref(j,:);
                            end
                            
                            u_Neumann = obj.M_DATA.bcNeu{k}(xlt,ylt,zlt,t,obj.M_DATA.param);
                            one       = ones(nof,nqn);
                            u_Neumann = u_Neumann.*one;
                            
                            x    =  obj.M_MESH.vertices(1,obj.M_MESH.boundaries(1:3, obj.M_MESH.Neumann_side{k}));
                            y    =  obj.M_MESH.vertices(2,obj.M_MESH.boundaries(1:3, obj.M_MESH.Neumann_side{k}));
                            z    =  obj.M_MESH.vertices(3,obj.M_MESH.boundaries(1:3, obj.M_MESH.Neumann_side{k}));
                            
                            areav = cross(  [x(2:3:end)-x(1:3:end);  y(2:3:end)-y(1:3:end);  z(2:3:end)-z(1:3:end)], ...
                                [x(3:3:end)-x(1:3:end);  y(3:3:end)-y(1:3:end);  z(3:3:end)-z(1:3:end)]);
                            
                            for l = 1 : nof
                                
                                area   = 0.5*norm(areav(:,l));
                                detjac = 2*area;
                                
                                face = obj.M_MESH.Neumann_side{k}(l);
                                
                                u_Neumann_loc  = u_Neumann(l,:).*wi;
                                u_Neumann_loc  = u_Neumann_loc(1,:)';
                                
                                Rrows(1+(l-1)*nbn:l*nbn)    = obj.M_MESH.boundaries(1:nbn,face);
                                Rcoef(1+(l-1)*nbn:l*nbn)    = detjac*phi*u_Neumann_loc;
                            end
                            F = F + sparse(Rrows+(k-1)*obj.M_MESH.numNodes,1,Rcoef,obj.M_MESH.dim*obj.M_MESH.numNodes,1);
                            
                        end
                    end
                    
                    % Pressure condition
                    for k = 1 : obj.M_MESH.dim
                        if ~isempty(obj.M_MESH.Pressure_side{k})
                            
                            [quad_points, wi] = quadrature(obj.M_MESH.dim-1, obj.M_FE_SPACE.quad_order);
                            csi = quad_points(1,:);
                            eta = quad_points(2,:);
                            [phi]          =  fem_basis(obj.M_MESH.dim, obj.M_FE_SPACE.fem, [csi; eta; 0*eta], 1);
                            eta1           =  1-csi-eta;
                            nqn            =  length(wi);
                            
                            nbn         = obj.M_MESH.numBoundaryDof;
                            
                            for flag = 1 : length(obj.M_DATA.flag_pressure{k})
                                
                                nof         = length(obj.M_MESH.Pressure_side_CompFlag{k,flag});
                                
                                Rrows       = zeros(nbn*nof,1);
                                Rcoef       = Rrows;
                                
                                xlt = zeros(nof,nqn); ylt = xlt; zlt = xlt;
                                coord_ref = [eta1; csi; eta];
                                for j = 1 : 2
                                    dof = obj.M_MESH.boundaries(j,obj.M_MESH.Pressure_side_CompFlag{k,flag});
                                    vtemp = obj.M_MESH.vertices(1,dof);
                                    xlt = xlt + vtemp'*coord_ref(j,:);
                                    vtemp = obj.M_MESH.vertices(2,dof);
                                    ylt = ylt + vtemp'*coord_ref(j,:);
                                    vtemp = obj.M_MESH.vertices(3,dof);
                                    zlt = zlt + vtemp'*coord_ref(j,:);
                                end
                                
                                if length(obj.M_DATA.bcPrex) == 1
                                    pressure = obj.M_DATA.bcPrex(xlt,ylt,zlt,t,obj.M_DATA.param);
                                else
                                    pressure = obj.M_DATA.bcPrex{obj.M_DATA.flag_pressure{k}(flag)}(xlt,ylt,zlt,t,param);
                                end
                                
                                one      = ones(nof,nqn);
                                pressure = pressure.*one;
                                
                                x    =  obj.M_MESH.vertices(1,obj.M_MESH.boundaries(1:3, obj.M_MESH.Pressure_side_CompFlag{k,flag}));
                                y    =  obj.M_MESH.vertices(2,obj.M_MESH.boundaries(1:3, obj.M_MESH.Pressure_side_CompFlag{k,flag}));
                                z    =  obj.M_MESH.vertices(3,obj.M_MESH.boundaries(1:3, obj.M_MESH.Pressure_side_CompFlag{k,flag}));
                                
                                areav = cross(  [x(2:3:end)-x(1:3:end);  y(2:3:end)-y(1:3:end);  z(2:3:end)-z(1:3:end)], ...
                                    [x(3:3:end)-x(1:3:end);  y(3:3:end)-y(1:3:end);  z(3:3:end)-z(1:3:end)]);
                                
                                for l = 1 : nof
                                    
                                    area   = 0.5*norm(areav(:,l));
                                    detjac = 2*area;
                                    
                                    face = obj.M_MESH.Pressure_side_CompFlag{k,flag}(l);
                                    
                                    pressure_loc  = pressure(l,:).*wi;
                                    pressure_loc  = pressure_loc(1,:)';
                                    
                                    Rrows(1+(l-1)*nbn:l*nbn)    = obj.M_MESH.boundaries(1:nbn,face);
                                    Rcoef(1+(l-1)*nbn:l*nbn)    = obj.M_MESH.Normal_Faces(k,face)*detjac*phi*pressure_loc;
                                end
                                F = F + sparse(Rrows+(k-1)*obj.M_MESH.numNodes,1,Rcoef,obj.M_MESH.dim*obj.M_MESH.numNodes,1);
                            end
                        end
                    end
                
            end
            
        end
        
        %==========================================================================
        %% Compute Surface Forces
        function F = compute_follower_Pload( obj, U_h, t )
            
            % only for Dim = 2
            
            if nargin < 3 || isempty(t)
                t = [];
            end
            
            F = sparse(obj.M_MESH.numNodes*obj.M_MESH.dim, 1);
           
            disp_nodes = [];
            for hh = 1 : obj.M_MESH.dim
               disp_nodes = [disp_nodes; U_h(1+(hh-1)*obj.M_FE_SPACE.numDofScalar:hh*obj.M_FE_SPACE.numDofScalar)']; 
            end
            def_nodes = obj.M_MESH.nodes + disp_nodes;
            def_vertices = def_nodes(:, 1:obj.M_MESH.numVertices);
            
            Normal_Faces =  ComputeSurfaceNormals3D(obj.M_MESH.boundaries(1:3,:),def_vertices, obj.M_MESH.elements(1:4,:));
            
            for k = 1 : obj.M_MESH.dim
                if ~isempty(obj.M_MESH.Pressure_side{k})
                    
                    [quad_points, wi] = quadrature(obj.M_MESH.dim-1, obj.M_FE_SPACE.quad_order);
                    csi = quad_points(1,:);
                    eta = quad_points(2,:);
                    [phi]          =  fem_basis(obj.M_MESH.dim, obj.M_FE_SPACE.fem, [csi; eta; 0*eta], 1);
                    eta1           =  1-csi-eta;
                    nqn            =  length(wi);
                    
                    nof         = length(obj.M_MESH.Pressure_side{k});
                    nbn         = obj.M_MESH.numBoundaryDof;
                    
                    Rrows       = zeros(nbn*nof,1);
                    Rcoef       = Rrows;
                    
                    xlt = zeros(nof,nqn); ylt = xlt; zlt = xlt;
                    coord_ref = [eta1; csi; eta];
                    for j = 1 : 2
                        dof = obj.M_MESH.boundaries(j,obj.M_MESH.Pressure_side{k});
                        vtemp = def_vertices(1,dof);
                        xlt = xlt + vtemp'*coord_ref(j,:);
                        vtemp = def_vertices(2,dof);
                        ylt = ylt + vtemp'*coord_ref(j,:);
                        vtemp = def_vertices(3,dof);
                        zlt = zlt + vtemp'*coord_ref(j,:);
                    end
                    
                    pressure = obj.M_DATA.bcPrex(xlt,ylt,zlt,t,obj.M_DATA.param);
                    one       = ones(nof,nqn);
                    pressure = pressure.*one;
                    
                    x    =  def_vertices(1,obj.M_MESH.boundaries(1:3, obj.M_MESH.Pressure_side{k}));
                    y    =  def_vertices(2,obj.M_MESH.boundaries(1:3, obj.M_MESH.Pressure_side{k}));
                    z    =  def_vertices(3,obj.M_MESH.boundaries(1:3, obj.M_MESH.Pressure_side{k}));
                    
                    areav = cross(  [x(2:3:end)-x(1:3:end);  y(2:3:end)-y(1:3:end);  z(2:3:end)-z(1:3:end)], ...
                        [x(3:3:end)-x(1:3:end);  y(3:3:end)-y(1:3:end);  z(3:3:end)-z(1:3:end)]);
                    
                    for l = 1 : nof
                        
                        area   = 0.5*norm(areav(:,l));
                        detjac = 2*area;
                        
                        face = obj.M_MESH.Pressure_side{k}(l);
                        
                        pressure_loc  = pressure(l,:).*wi;
                        pressure_loc  = pressure_loc(1,:)';
                        
                        Rrows(1+(l-1)*nbn:l*nbn)    = obj.M_MESH.boundaries(1:nbn,face);
                        Rcoef(1+(l-1)*nbn:l*nbn)    = Normal_Faces(k,face)*detjac*phi*pressure_loc;
                        
                    end
                    F = F + sparse(Rrows+(k-1)*obj.M_MESH.numNodes,1,Rcoef,obj.M_MESH.dim*obj.M_MESH.numNodes,1);
                    
                end
            end
            
        end
        
        %==========================================================================
        %% Compute External Forces
        function F_ext = compute_external_forces( obj, t )
            
            if nargin < 2 || isempty(t)
                t = [];
            end
            
            F_ext = compute_volumetric_forces( obj, t ) + compute_surface_forces( obj, t );
            
        end
        
        %==========================================================================
        %% Compute mass matrix
        function [M] = compute_mass( obj )
            
            % C_OMP assembly, returns matrices in sparse vector format
            [rowM, colM, coefM] = Mass_assembler_C_omp(obj.M_MESH.dim, obj.M_MESH.elements, obj.M_FE_SPACE.numElemDof, ...
                obj.M_FE_SPACE.quad_weights, obj.M_MESH.jac, obj.M_FE_SPACE.phi);
            
            % Build sparse matrix
            M_scalar   = GlobalAssemble(rowM, colM, coefM, obj.M_MESH.numNodes, obj.M_MESH.numNodes);
            M          = [];
            for k = 1 : obj.M_FE_SPACE.numComponents
                M = blkdiag(M, M_scalar);
            end
            
        end
        
        %==========================================================================
        %% Compute stress
        function [P, Sigma] = compute_stress(obj, U_h)
            % C_OMP compute element stresses, returns dense matrix of size
            % N_elements x MESH.dim^2
            %
            % P is the First Piola-Kirchhof stress tensor, while Sigma is
            % the Cauchy stress tensor ( Sigma = 1 / det(F) * P * F^T )
            
            [quad_nodes, quad_weights]   = quadrature(obj.M_MESH.dim, 1);
            [phi, dphi_ref]              = fem_basis(obj.M_FE_SPACE.dim, obj.M_FE_SPACE.fem, quad_nodes);
            
            [P, Sigma] = CSM_assembler_C_omp(obj.M_MESH.dim, [obj.M_MaterialModel,'_stress'], obj.M_MaterialParam, ...
                U_h, obj.M_MESH.elements, obj.M_FE_SPACE.numElemDof,...
                quad_weights, obj.M_MESH.invjac, obj.M_MESH.jac, phi, dphi_ref);
            
        end
        
        %==========================================================================
        %% Compute internal forces
        function [F_in] = compute_internal_forces(obj, U_h)
            
            % C_OMP assembly, returns matrices in sparse vector format
            [rowG, coefG] = ...
                CSM_assembler_C_omp(obj.M_MESH.dim, [obj.M_MaterialModel,'_forces'], obj.M_MaterialParam, full( U_h ), ...
                obj.M_MESH.elements, obj.M_FE_SPACE.numElemDof, ...
                obj.M_FE_SPACE.quad_weights, obj.M_MESH.invjac, obj.M_MESH.jac, obj.M_FE_SPACE.phi, obj.M_FE_SPACE.dphi_ref);
            
            % Build sparse matrix and vector
            F_in    = GlobalAssemble(rowG, 1, coefG, obj.M_MESH.numNodes*obj.M_MESH.dim, 1);
            
        end
        
        %==========================================================================
        %% Compute internal forces Jacobian
        function [dF_in] = compute_jacobian(obj, U_h)
            
            if nargin < 2 || isempty(U_h)
                U_h = zeros(obj.M_MESH.dim*obj.M_MESH.numNodes,1);
            end

            % C_OMP assembly, returns matrices in sparse vector format
            [rowdG, coldG, coefdG] = ...
                CSM_assembler_C_omp(obj.M_MESH.dim, [obj.M_MaterialModel,'_jacobian'], obj.M_MaterialParam, full( U_h ), ...
                obj.M_MESH.elements, obj.M_FE_SPACE.numElemDof, ...
                obj.M_FE_SPACE.quad_weights, obj.M_MESH.invjac, obj.M_MESH.jac, obj.M_FE_SPACE.phi, obj.M_FE_SPACE.dphi_ref);
            
            % Build sparse matrix and vector
            dF_in   = GlobalAssemble(rowdG, coldG, coefdG, obj.M_MESH.numNodes*obj.M_MESH.dim, obj.M_MESH.numNodes*obj.M_MESH.dim);
        end
        
        %==========================================================================
        %% Compute Prestress vector and jacobian
        function [R_P, J_P, S_np1] = compute_prestress(obj, U_h, S_0)
            
            if nargin < 2 || isempty(U_h)
                U_h = zeros(obj.M_MESH.dim*obj.M_MESH.numNodes,1);
            end
            
            if nargin < 3 || isempty(S_0)
                S_0 = zeros(obj.M_MESH.numElem*length(obj.M_FE_SPACE.quad_weights)*obj.M_MESH.dim*obj.M_MESH.dim, 1);
            end

            % C_OMP assembly, returns matrices in sparse vector format
            [rowdG, coldG, coefdG, rowG, coefG, S_np1] = ...
                CSM_assembler_C_omp(obj.M_MESH.dim, [obj.M_MaterialModel,'_prestress'], obj.M_MaterialParam, full( U_h ), ...
                obj.M_MESH.elements, obj.M_FE_SPACE.numElemDof, ...
                obj.M_FE_SPACE.quad_weights, obj.M_MESH.invjac, obj.M_MESH.jac, obj.M_FE_SPACE.phi, obj.M_FE_SPACE.dphi_ref, S_0);
            
            % Build sparse matrix and vector
            R_P   = GlobalAssemble(rowG, 1, coefG, obj.M_MESH.numNodes*obj.M_MESH.dim, 1);
            J_P   = GlobalAssemble(rowdG, coldG, coefdG, obj.M_MESH.numNodes*obj.M_MESH.dim, obj.M_MESH.numNodes*obj.M_MESH.dim);
        end
        
        %==========================================================================
        %% Assemble Robin Condition: Pn + K d = 0 on \Gamma_Robin, with K = ElasticCoefRobin
        function [A] = assemble_ElasticRobinBC(obj)
            
            A = sparse(obj.M_MESH.numNodes*obj.M_MESH.dim, obj.M_MESH.numNodes*obj.M_MESH.dim);
            
            for k = 1 : obj.M_MESH.dim
                if ~isempty(obj.M_MESH.Robin_side{k})
                    
                    [quad_points, wi] = quadrature(obj.M_MESH.dim-1, obj.M_FE_SPACE.quad_order);
                    csi = quad_points(1,:);
                    eta = quad_points(2,:);
                    [phi]       =  fem_basis(obj.M_MESH.dim, obj.M_FE_SPACE.fem, [csi; eta; 0*eta], 1);
                    
                    nof         = length(obj.M_MESH.Robin_side{k});
                    nbn         = obj.M_MESH.numBoundaryDof;
                    
                    Arows       = zeros(nbn*nbn*nof,1);
                    Acols       = Arows;
                    Acoef       = Arows;
                    [rows,cols] = meshgrid(1:nbn,1:nbn);
                    rows        = rows(:);
                    cols        = cols(:);
                    
                    x    =  obj.M_MESH.vertices(1,obj.M_MESH.boundaries(1:3, obj.M_MESH.Robin_side{k}));
                    y    =  obj.M_MESH.vertices(2,obj.M_MESH.boundaries(1:3, obj.M_MESH.Robin_side{k}));
                    z    =  obj.M_MESH.vertices(3,obj.M_MESH.boundaries(1:3, obj.M_MESH.Robin_side{k}));
                    
                    areav = cross(  [x(2:3:end)-x(1:3:end);  y(2:3:end)-y(1:3:end);  z(2:3:end)-z(1:3:end)], ...
                        [x(3:3:end)-x(1:3:end);  y(3:3:end)-y(1:3:end);  z(3:3:end)-z(1:3:end)]);
                    
                    MASS_loc = (phi.*repmat(wi,nbn,1))*phi';
                    MASS_loc = MASS_loc(:);
                    for l = 1 : nof
                        
                        area   = 0.5*norm(areav(:,l));
                        detjac = 2*area;
                        
                        face = obj.M_MESH.Robin_side{k}(l);
                        
                        Arows(1+(l-1)*nbn*nbn:l*nbn*nbn)   =  obj.M_MESH.boundaries(rows,face);
                        Acols(1+(l-1)*nbn*nbn:l*nbn*nbn)   =  obj.M_MESH.boundaries(cols,face);
                        Acoef(1+(l-1)*nbn*nbn:l*nbn*nbn)   =  detjac*obj.M_DATA.ElasticCoefRobin*MASS_loc;
                        
                    end
                    A = A + sparse(Arows+(k-1)*obj.M_MESH.numNodes,Acols+(k-1)*obj.M_MESH.numNodes,...
                        Acoef,obj.M_MESH.dim*obj.M_MESH.numNodes,obj.M_MESH.dim*obj.M_MESH.numNodes);
                    
                end
            end
        end
        
    end
    
end