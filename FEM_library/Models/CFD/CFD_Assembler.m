function [varargout] = CFD_Assembler(output, MESH, DATA, FE_SPACE_v, FE_SPACE_p, U_h, t, subdomain)
%CFD_ASSEMBLER assembler for 2D/3D CFD models

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri@epfl.ch>

if nargin < 6 || isempty(U_h)
    U_h = zeros(FE_SPACE_v.numDof + FE_SPACE_p.numDof,1);
end

if nargin < 7
    t = [];
end

if nargin < 8
    subdomain = [];
end

if ~isempty(subdomain)
    index_subd = [];
    for q = 1 : length(subdomain)
        index_subd = [index_subd find(MESH.elements(FE_SPACE.numElemDof+1,:) == subdomain(q))];
    end
    MESH.elements = MESH.elements(:,index_subd);
    MESH.numElem  = size(MESH.elements,2);
else
    index_subd = [1:MESH.numElem];
end

switch output
    
    case 'volume_force'
        
        F_volume = compute_external_forces(MESH, DATA, FE_SPACE_v, t, index_subd);
        
        varargout{1} = [F_volume; zeros(FE_SPACE_p.numDof,1)];
        
    case 'Stokes'
        
        [A] = compute_Stokes_matrix(DATA.kinematic_viscosity, MESH, FE_SPACE_v, FE_SPACE_p, index_subd);
        
        varargout{1} = A;    
        
    case 'convective_Oseen'
        
        [C] = compute_convective_Oseen_matrix(DATA.density, MESH, FE_SPACE_v, FE_SPACE_p, U_h, index_subd);
                
        varargout{1} = 0*C;
        
    case 'convective'
        
        [C1, C2] = compute_convective_matrix(DATA.density, MESH, FE_SPACE_v, FE_SPACE_p, U_h, index_subd);
        
        varargout{1} = C1;
        varargout{2} = C2;
        
        
    case 'mass_velocity'
        
        [M_v] = compute_mass( MESH, FE_SPACE_v, index_subd);
        
        varargout{1} = M_v;
        
    case 'mass_pressure'
        
        [M_p] = compute_mass( MESH, FE_SPACE_p, index_subd);
        
        varargout{1} = M_p;

        
    otherwise
        error('output option not available')
end


end

%==========================================================================
function [F_ext] = compute_external_forces(MESH, DATA, FE_SPACE, t, index_subd)

% Computations of all quadrature nodes in the elements
coord_ref = MESH.chi;
switch MESH.dim
    
    case 2
        
        x = zeros(MESH.numElem,FE_SPACE.numQuadNodes); y = x;
        for j = 1 : 3
            i = MESH.elements(j,:);
            vtemp = MESH.vertices(1,i);
            x = x + vtemp'*coord_ref(j,:);
            vtemp = MESH.vertices(2,i);
            y = y + vtemp'*coord_ref(j,:);
        end
        
        % Evaluation of external forces in the quadrature nodes
        for k = 1 : MESH.dim
            f{k}  = DATA.force{k}(x,y,t,DATA.param);
        end
        
    case 3
        
        x = zeros(MESH.numElem,FE_SPACE.numQuadNodes); y = x; z = x;
        
        for j = 1 : 4
            i = MESH.elements(j,:);
            vtemp = MESH.vertices(1,i);
            x = x + vtemp'*coord_ref(j,:);
            vtemp = MESH.vertices(2,i);
            y = y + vtemp'*coord_ref(j,:);
            vtemp = MESH.vertices(3,i);
            z = z + vtemp'*coord_ref(j,:);
        end
        
        % Evaluation of external forces in the quadrature nodes
        for k = 1 : MESH.dim
            f{k}  = DATA.force{k}(x,y,z,t,DATA.param);
        end
        
end
% C_OMP assembly, returns matrices in sparse vector format

F_ext = [];
for k = 1 : MESH.dim
    
    [rowF, coefF] = CSM_assembler_ExtForces(f{k}, MESH.elements, FE_SPACE.numElemDof, ...
        FE_SPACE.quad_weights, MESH.jac(index_subd), FE_SPACE.phi);
    
    % Build sparse matrix and vector
    F_ext    = [F_ext; sparse(rowF, 1, coefF, FE_SPACE.numDofScalar, 1)];
    
end

end
%==========================================================================
function [M] = compute_mass( MESH, FE_SPACE, index_subd)

% C_OMP assembly, returns matrices in sparse vector format
[rowM, colM, coefM] = Mass_assembler_C_omp(MESH.dim, MESH.elements, FE_SPACE.numElemDof, ...
     FE_SPACE.quad_weights, MESH.jac(index_subd), FE_SPACE.phi);

% Build sparse matrix
M_scalar   = sparse(rowM, colM, coefM, FE_SPACE.numDofScalar, FE_SPACE.numDofScalar);
M          = [];
for k = 1 : FE_SPACE.numComponents
    M = blkdiag(M, M_scalar);
end

end
%==========================================================================
function [A] =  compute_Stokes_matrix(viscosity, MESH, FE_SPACE_v, FE_SPACE_p, index_subd)
        
% C_OMP assembly, returns matrices in sparse vector format
[rowA, colA, coefA] = ...
    CFD_assembler_C_omp('Stokes', viscosity, MESH.dim, MESH.elements, ...
    FE_SPACE_v.numElemDof, FE_SPACE_p.numElemDof, FE_SPACE_v.numDof, ...
    FE_SPACE_v.quad_weights, MESH.invjac(index_subd,:,:), MESH.jac(index_subd), ...
    FE_SPACE_v.phi, FE_SPACE_v.dphi_ref, FE_SPACE_p.phi);

% Build sparse matrix
A   = sparse(rowA, colA, coefA, FE_SPACE_v.numDof + FE_SPACE_p.numDof, FE_SPACE_v.numDof + FE_SPACE_p.numDof);

end
%==========================================================================
function [C] =  compute_convective_Oseen_matrix(density, MESH, FE_SPACE_v, FE_SPACE_p, U_h, index_subd)
        
% C_OMP assembly, returns matrices in sparse vector format
[rowA, colA, coefA] = ...
    CFD_assembler_C_omp('convective_Oseen', density, MESH.dim, MESH.elements, ...
    FE_SPACE_v.numElemDof, FE_SPACE_v.numDof, ...
    FE_SPACE_v.quad_weights, MESH.invjac(index_subd,:,:), MESH.jac(index_subd), ...
    FE_SPACE_v.phi, FE_SPACE_v.dphi_ref, U_h);

% Build sparse matrix

C   = sparse(rowA, colA, coefA, FE_SPACE_v.numDof + FE_SPACE_p.numDof, FE_SPACE_v.numDof + FE_SPACE_p.numDof);
C   = density * C;
end
%==========================================================================
function [C1, C2] =  compute_convective_matrix(density, MESH, FE_SPACE_v, FE_SPACE_p, U_h, index_subd)
        
% C_OMP assembly, returns matrices in sparse vector format
[rowA, colA, coefA, rowB, colB, coefB] = ...
    CFD_assembler_C_omp('convective', density, MESH.dim, MESH.elements, ...
    FE_SPACE_v.numElemDof, FE_SPACE_v.numDof, ...
    FE_SPACE_v.quad_weights, MESH.invjac(index_subd,:,:), MESH.jac(index_subd), ...
    FE_SPACE_v.phi, FE_SPACE_v.dphi_ref, U_h);

% Build sparse matrix

C1   = sparse(rowA, colA, coefA, FE_SPACE_v.numDof + FE_SPACE_p.numDof, FE_SPACE_v.numDof + FE_SPACE_p.numDof);
C2   = sparse(rowB, colB, coefB, FE_SPACE_v.numDof + FE_SPACE_p.numDof, FE_SPACE_v.numDof + FE_SPACE_p.numDof);

end
%==========================================================================