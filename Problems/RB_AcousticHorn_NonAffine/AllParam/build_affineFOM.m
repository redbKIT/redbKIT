function [ FOM ] = build_affineFOM(elements, vertices, boundaries, fem, data_file, tolPOD)
%BUILD_AFFINEFOM returns a FOM (for full-order model) struct containing all
%the structures required for the reduction.

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

FOM.model = 'ADR';

%% Parameters
FOM.P  = 5;
% Full variations
FOM.mu_min = [ 50    -0.03     -0.03     -0.03    -0.03];
FOM.mu_max = [ 1000  0.03       0.03      0.03     0.03];
FOM.mu_bar = [ 500    0          0         0        0]; 

%% Read problem parameters and BCs from data_file
DATA       = read_DataFile(data_file);
DATA.param = zeros(1,5);
DATA.shape_param = [2:5];
t          = [];

%% Set quad_order
quad_order                  = 4;  

%% Fill MESH data structure
dim              = 2;
[ MESH ] = buildMESH( dim, elements, vertices, boundaries, fem, quad_order, DATA );

%% Create and fill the FE_SPACE data structure
[ FE_SPACE ] = buildFESpace( MESH, fem, 1, quad_order );

fprintf('\n **** PROBLEM''S SIZE INFO ****\n');
fprintf(' * Number of Vertices  = %d \n',MESH.numVertices);
fprintf(' * Number of Elements  = %d \n',MESH.numElem);
fprintf(' * Number of Nodes     = %d \n',MESH.numNodes);
fprintf(' * Number of Dofs      = %d \n',FE_SPACE.numDof);
fprintf('-------------------------------------------\n');

%% Fill FOM structures
FOM.MESH     = MESH;
FOM.FE_SPACE = FE_SPACE;
FOM.DATA     = DATA;

%% Assemble MATRIX wih DEIM
mu_train_Dimension = 250; % number of samples
mu_cube            = lhsdesign(mu_train_Dimension,FOM.P); % generate normalized design
mu_train           = bsxfun(@plus,FOM.mu_min,bsxfun(@times,mu_cube,(FOM.mu_max-FOM.mu_min)));


%% compute FEM snapshots, store matrices and rhs
ndf       = FE_SPACE.numDof; % number of DOFs
S_matrix  = sparse(ndf*ndf,mu_train_Dimension);
S_rhs     = sparse(ndf,mu_train_Dimension);

parfor i = 1 : size(mu_train,1)
    
    MESH_tmp            =  MESH;
    MESH_tmp.vertices   =  RBF_MeshDeformation(mu_train(i,FOM.DATA.shape_param),MESH.vertices);
    [MESH_tmp.jac, MESH_tmp.invjac, MESH_tmp.h] = geotrasf(MESH_tmp.dim, MESH_tmp.vertices, MESH_tmp.elements);   
    DATA_tmp            =  DATA;
    DATA_tmp.param      =  mu_train(i,:);
    
    [A, F]              =  ADR_Assembler(MESH_tmp, DATA_tmp, FE_SPACE);
    [A_in, F_in]        =  ADR_ApplyBC(A, F, FE_SPACE, MESH_tmp, DATA_tmp);

    S_matrix(:,i)       = A_in(:);
    S_rhs(:,i)          = F_in;
    
    fprintf('\nAssembled Snapshot Matrix and Vector %d of %d', i, mu_train_Dimension);
end

%% Generate Approximate Affine Model by (M)DEIM
[FOM.HyRED, U_matrix, U_rhs]  = ADR_HyperReduction_offline(FOM, S_rhs, S_matrix, tolPOD);
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
DATA.diffusion =  @(x,y,t,param)(1+0.*x.*y);
X              =  ADR_Assembler(MESH, DATA, FE_SPACE, 'diffusion', [], [], []);
DATA.reaction  =  @(x,y,t,param)(1+0.*x.*y);
M              =  ADR_Assembler(MESH, DATA, FE_SPACE, 'reaction');

FOM.Xnorm =  X(FOM.MESH.internal_dof,FOM.MESH.internal_dof) +...
             M(FOM.MESH.internal_dof,FOM.MESH.internal_dof);

end
