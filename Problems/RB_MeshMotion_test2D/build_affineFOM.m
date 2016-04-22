function [ FOM ] = build_affineFOM(elements, vertices, boundaries, fem, data_file, tolPOD)
%BUILD_AFFINEFOM returns a FOM (for full-order model) struct containing all
%the structures required for the reduction.

%   Author: F. Negri (federico.negri@epfl.ch) 2015
%   Copyright (C) Federico Negri, CMCS, EPFL

dim              = 2;

FOM.model   = 'CSM';
%% Parameters
FOM.P  = 3;

% Only rhs
FOM.mu_min = [ -0.25   -0.25  0.0   ];
FOM.mu_max = [  0.25    0.25  1.2   ];
FOM.mu_bar = [   0     0  0.0   ]; 

%% Read problem parameters and BCs from data_file
DATA       = CSM_read_DataFile(data_file, dim, [0  0 0]);
DATA.param = zeros(1,3);
t          = [];

%% Set quad_order
quad_order                  = 4;  

%% Fill MESH data structure
[ MESH ] = buildMESH( dim, elements, vertices, boundaries, fem, quad_order, DATA, 'CSM' );

%% Create and fill the FE_SPACE data structure
[ FE_SPACE ] = buildFESpace( MESH, fem, 1, quad_order );

fprintf('\n **** PROBLEM''S SIZE INFO ****\n');
fprintf(' * Number of Vertices  = %d \n',MESH.numVertices);
fprintf(' * Number of Elements  = %d \n',MESH.numElem);
fprintf(' * Number of Nodes     = %d \n',MESH.numNodes);
fprintf(' * Number of Dofs      = %d \n',length(MESH.internal_dof));
fprintf('-------------------------------------------\n');

%% Fill FOM structures
FOM.MESH     = MESH;
FOM.FE_SPACE = FE_SPACE;
FOM.DATA     = DATA;

%% Assemble MATRIX wih DEIM
mu_train_Dimension = 30; % number of samples
mu_cube            = lhsdesign(mu_train_Dimension,FOM.P); % generate normalized design
mu_train           = bsxfun(@plus,FOM.mu_min,bsxfun(@times,mu_cube,(FOM.mu_max-FOM.mu_min)));


%% compute FEM snapshots, store matrices and rhs
ndf       = FE_SPACE.numDof; % number of DOFs
S_matrix  = sparse(ndf*ndf,mu_train_Dimension);
S_rhs     = sparse(ndf,mu_train_Dimension);

parfor i = 1 : size(mu_train,1)
    
    DATA_tmp            =  CSM_read_DataFile(data_file, dim, mu_train(i,:));
    DATA_tmp.param      =  mu_train(i,:);
    
    [~, A]              =  CSM_Assembler('internal_forces', MESH, DATA_tmp, FE_SPACE);
    [A_in, F_in]        =  CSM_ApplyBC(A, [], FE_SPACE, MESH, DATA_tmp);
    S_matrix(:,i)       = A_in(:);
    S_rhs(:,i)          = F_in;
    
    fprintf('\nAssembled Snapshot Matrix and Vector %d of %d', i, mu_train_Dimension);
end

%% Generate Approximate Affine Model by (M)DEIM
[FOM.HyRED, U_matrix, U_rhs]  = CSM_HyperReduction_offline(FOM, S_rhs, S_matrix, tolPOD);
clear S_matrix S_rhs;

FOM.Qa    = size(U_matrix,2);
FOM.Qf    = size(U_rhs,2);

for q = 1 : FOM.Qa
      FOM.Aq{q} = reshape(U_matrix(:,q), ndf, ndf);
end

for q = 1 : FOM.Qf
      FOM.Fq{q} = U_rhs(:,q);
end

%% Compute Xnorm
DATA.diffusion    =  @(x,y,t,param)(1+0.*x.*y);
DATA.transport{1} =  @(x,y,t,param)(0+0.*x.*y);
DATA.transport{2} =  @(x,y,t,param)(0+0.*x.*y);
DATA.force        =  @(x,y,t,param)(0+0.*x.*y);
DATA.reaction  =  @(x,y,t,param)(1+0.*x.*y);
X              =  ADR_Assembler(MESH, DATA, FE_SPACE, 'diffusion', [], [], []);

X = blkdiag(X,X);
FOM.Xnorm =  X(MESH.internal_dof, MESH.internal_dof);

end
