function [A, F, M] = ADR_Assembler(MESH, DATA, FE_SPACE, OPERATOR, TC_d, TC_t, subdomain, t, stabilization, dt)
%ADR_ASSEMBLER assembler for 2D/3D ADR equations
%
%   [A, F] = ADR_ASSEMBLER(MESH, DATA, FE_SPACE) given MESH and 
%   FE_SPACE structs, assemble the stiffness matrix A and rhs vector F
%   of the advection-diffusion-reaction differential operator using the
%   parameters specified in the DATA struct. The resulting matrix and
%   vector doesn't take into account for Boundary Conditions.
%
%   [A, F, M] = ADR_ASSEMBLER(MESH, DATA, FE_SPACE, T) as before, but also
%   returns the mass matrix M. Moreover, the problem's parameters are
%   evaluated with respect to a given time T.
%
%   [A] = ADR_ASSEMBLER(MESH, DATA, FE_SPACE, 'diffusion', [i j], [], k);
%   returns the matrix corresponding to the (i,j)-th second derivative
%   (diffusion operator) assembled over the subdomain \Omega_k
%   A spatial dependent coefficient has to be defined in DATA.diffusion
%
%   [A] = ADR_ASSEMBLER(MESH, DATA, FE_SPACE, 'transport', [], [i], k);
%   returns the matrix corresponding to the (i)-th first derivative 
%   (advection operator) assembled over the subdomain \Omega_k
%   A spatial dependent coefficient has to be defined in DATA.transport{i}

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri@epfl.ch> 

if nargin < 4 || isempty(OPERATOR)
    OPERATOR = 'all';
end

if nargin < 5 || isempty(TC_d)
    TC_d = [10 10];% diagonal components of the diffusion operator by default
end

if nargin < 6 || isempty(TC_t)
    TC_t = 10; % all components of the transport operator by default
end

if nargin < 7
    subdomain = [];
end

if nargin < 8
    t = [];
end

if nargin < 9
    stabilization = [];
end

if nargin < 10
    dt = 1;
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


%% Computations of all quadrature nodes in the elements
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
        
        %% Evaluation of coefficients in the quadrature nodes
        mu  = EvalDataQuad(MESH, DATA.param, FE_SPACE, t, {x, y}, DATA.diffusion, index_subd);
        si  = EvalDataQuad(MESH, DATA.param, FE_SPACE, t, {x, y}, DATA.reaction, index_subd);
        f   = EvalDataQuad(MESH, DATA.param, FE_SPACE, t, {x, y}, DATA.force, index_subd);
        
        bx  = EvalDataQuad(MESH, DATA.param, FE_SPACE, t, {x, y}, DATA.transport{1}, index_subd);
        by  = EvalDataQuad(MESH, DATA.param, FE_SPACE, t, {x, y}, DATA.transport{2}, index_subd);
        
        b = [bx by];
        
    case 3
        x = zeros(MESH.numElem,FE_SPACE.numQuadNodes); y = x; z = x;

        for j = 1:4
            i = MESH.elements(j,:);
            vtemp = MESH.vertices(1,i);
            x = x + vtemp'*coord_ref(j,:);
            vtemp = MESH.vertices(2,i);
            y = y + vtemp'*coord_ref(j,:);
            vtemp = MESH.vertices(3,i);
            z = z + vtemp'*coord_ref(j,:);
        end
        
        %% Evaluation of coefficients in the quadrature nodes
        mu  = EvalDataQuad(MESH, DATA.param, FE_SPACE, t, {x, y, z}, DATA.diffusion, index_subd);
        si  = EvalDataQuad(MESH, DATA.param, FE_SPACE, t, {x, y, z}, DATA.reaction, index_subd);
        f   = EvalDataQuad(MESH, DATA.param, FE_SPACE, t, {x, y, z}, DATA.force, index_subd);
  
        bx  = EvalDataQuad(MESH, DATA.param, FE_SPACE, t, {x, y, z}, DATA.transport{1}, index_subd);
        by  = EvalDataQuad(MESH, DATA.param, FE_SPACE, t, {x, y, z}, DATA.transport{2}, index_subd);
        bz  = EvalDataQuad(MESH, DATA.param, FE_SPACE, t, {x, y, z}, DATA.transport{3}, index_subd);
        
        b = [bx by bz];
end

%% Assembly
if isempty( stabilization )
    
    % C assembly, returns matrices in sparse vector format
    [Arows, Acols, Acoef, Mcoef, Rrows, Rcoef] = ADR_assembler_C_omp(MESH.dim, OPERATOR, TC_d, TC_t, MESH.elements, FE_SPACE.numElemDof, mu, b, si, f,...
        FE_SPACE.quad_weights, MESH.invjac(index_subd,:,:), MESH.jac(index_subd), FE_SPACE.phi, FE_SPACE.dphi_ref);
    
    % Build sparse matrices and rhs
    A    = GlobalAssemble(Arows,Acols,Acoef,MESH.numNodes,MESH.numNodes);
    M    = GlobalAssemble(Arows,Acols,Mcoef,MESH.numNodes,MESH.numNodes);
    F    = GlobalAssemble(Rrows,1,Rcoef,MESH.numNodes,1);
    
else
        
    switch stabilization
       
        case 'SUPG'
            
            % C assembly, returns matrices in sparse vector format
            [Arows, Acols, Acoef, Mcoef, Rrows, Rcoef] = ADR_SUPGassembler_C_omp(MESH.dim, stabilization, dt, MESH.elements, FE_SPACE.numElemDof, mu, b, si, f,...
                FE_SPACE.quad_weights, MESH.invjac(index_subd,:,:), MESH.jac(index_subd), FE_SPACE.phi, FE_SPACE.dphi_ref);
            
            M = [];
            
        case 'SUPGt'
            
            % C assembly, returns matrices in sparse vector format
            [Arows, Acols, Acoef, Mcoef, Rrows, Rcoef] = ADR_SUPGassembler_C_omp(MESH.dim, stabilization, dt, MESH.elements, FE_SPACE.numElemDof, mu, b, si, f,...
                FE_SPACE.quad_weights, MESH.invjac(index_subd,:,:), MESH.jac(index_subd), FE_SPACE.phi, FE_SPACE.dphi_ref);
        
            M    = GlobalAssemble(Arows,Acols,Mcoef,MESH.numNodes,MESH.numNodes);
    end
    
    % Build sparse matrices and rhs
    A    = GlobalAssemble(Arows,Acols,Acoef,MESH.numNodes,MESH.numNodes);
    F    = GlobalAssemble(Rrows,1,Rcoef,MESH.numNodes,1);
    
end

return