function [ FOM ] = build_affineFOM(elements, vertices, boundaries, fem, data_file, tolPOD)
%BUILD_AFFINEFOM returns a FOM (for full-order model) struct containing all
%the structures required for the reduction.

%   Author: F. Negri (federico.negri@epfl.ch) 2015
%   Copyright (C) Federico Negri, CMCS, EPFL

%% Parameters
FOM.P  = 5;
FOM.mu_min = [ 50    -0.03     -0.03     -0.03    -0.03];
FOM.mu_max = [ 1000  0.03       0.03      0.03     0.03];
FOM.mu_bar = [ 500    0          0         0        0]; 

%% Read problem parameters and BCs from data_file
DATA       = read_DataFile(data_file);
DATA.param = zeros(1,5);
DATA.shape_param = [2:5];
t          = [];

%% Fill MESH data structure
MESH.dim         = 2;
MESH.vertices    = vertices;
MESH.boundaries  = boundaries;
MESH.elements    = elements;
MESH.numVertices = size(vertices,2);

%% Build higher order (P2 or P3) mesh if required
if ~strcmp(fem,'P1')
    [MESH.elements, MESH.nodes, MESH.boundaries] = ...
        feval(strcat('P1to',fem,'mesh','2D'),elements,vertices, boundaries);
else
    MESH.nodes = vertices;
end

%% Update Mesh data with BC information and geometrical maps
[numElemDof,numBoundaryDof]  = select2D(fem);
MESH.numNodes                = size(MESH.nodes,2);
MESH.numElem                 = size(MESH.elements,2);
MESH.numBoundaryDof          = numBoundaryDof;

% Update MESH with BC information
[MESH]         = BC_info(MESH, DATA);

% Compute geometrical map (ref to physical elements) information
[MESH.jac, MESH.invjac, MESH.h] = geotrasf2D(MESH.vertices, MESH.elements);   

% Compute quadrature nodes and weights on the reference element
quad_order                  = 4;  
[quad_nodes, quad_weights]  = dunavant_quad(quad_order);

% Evaluate P1 geometrical mapping basis functions in the quad points
[MESH.chi]                  =  fem_basis2D('P1', quad_nodes(1,:), quad_nodes(2,:));

%% Create and fill the FE_SPACE data structure
FE_SPACE.fem              = fem;
FE_SPACE.numDof           = length(MESH.internal_dof);
FE_SPACE.numElemDof       = numElemDof;
FE_SPACE.numBoundaryDof   = numBoundaryDof;

% Store quadrature nodes and weights on the reference element
FE_SPACE.quad_order    = quad_order;
FE_SPACE.quad_nodes    = quad_nodes;
FE_SPACE.quad_weights  = quad_weights;
FE_SPACE.numQuadNodes  = length(FE_SPACE.quad_nodes);

% Evaluate basis functions in the quad points on the reference element
[FE_SPACE.phi, FE_SPACE.dcsiphi, FE_SPACE.detaphi]  =  ...
    fem_basis2D(FE_SPACE.fem, FE_SPACE.quad_nodes(1,:), FE_SPACE.quad_nodes(2,:));


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
% [mu_train]         = FullFactorial_ParameterSpace(FOM.P, FOM.mu_min, FOM.mu_max, FOM.mu_ref, Tgrid_dimension);
% mu_train           = [mu_ref; mu_train];
% mu_train_Dimension = size(mu_train,1);

mu_train_Dimension = 100; % number of samples
mu_cube            = lhsdesign(mu_train_Dimension,FOM.P); % generate normalized design
mu_train           = bsxfun(@plus,FOM.mu_min,bsxfun(@times,mu_cube,(FOM.mu_max-FOM.mu_min)));


%% compute FEM snapshots, store matrices and rhs
ndf       = FE_SPACE.numDof; % number of DOFs
S_matrix  = sparse(ndf*ndf,mu_train_Dimension);
S_rhs     = sparse(ndf,mu_train_Dimension);

parfor i = 1 : size(mu_train,1)
    
    MESH_tmp            =  MESH;
    MESH_tmp.vertices   =  RBF_DeformGeometry(mu_train(i,FOM.DATA.shape_param),MESH.vertices);
    [MESH_tmp.jac, MESH_tmp.invjac, MESH_tmp.h] = geotrasf2D(MESH_tmp.vertices, MESH_tmp.elements);   
    DATA_tmp            =  DATA;
    DATA_tmp.param      =  mu_train(i,:);
    
    [A, F]              =  Assembler_2D(MESH_tmp, DATA_tmp, FE_SPACE);
    [A_in, F_in]        =  ApplyBC_2D(A, F, FE_SPACE, MESH_tmp, DATA_tmp);

    S_matrix(:,i)       = A_in(:);
    S_rhs(:,i)          = F_in;
    
    fprintf('\nAssembled Snapshot Matrix and Vector %d of %d', i, mu_train_Dimension);
end

%% Generate Approximate Affine Model by (M)DEIM
[FOM.HyRED, U_matrix, U_rhs]  = HyperReduction_offline(FOM, S_rhs, S_matrix, tolPOD);

FOM.Qa    = size(U_matrix,2);
FOM.Qf    = size(U_rhs,2);

for q = 1 : FOM.Qa
      FOM.Aq{q} = reshape(U_matrix(:,q), ndf, ndf);
end

for q = 1 : FOM.Qf
      FOM.Fq{q} = U_rhs(:,q);
end
%clear U_matrix U_rhs S_matrix S_rhs;

%% Compute Xnorm
DATA.diffusion =  @(x,y,t,param)(1+0.*x.*y);
X              =  Assembler_2D(MESH, DATA, FE_SPACE, 'diffusion', [], [], []);
DATA.reaction  =  @(x,y,t,param)(1+0.*x.*y);
M              =  Assembler_2D(MESH, DATA, FE_SPACE, 'reaction');

FOM.Xnorm =  X(FOM.MESH.internal_dof,FOM.MESH.internal_dof) +...
             M(FOM.MESH.internal_dof,FOM.MESH.internal_dof);

%clear MESH DATA FE_SPACE elements vertices boundaries;       
%FOM.u_D = @(x,mu) [];

end
