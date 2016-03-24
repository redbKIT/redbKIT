clear all
clc
close all

addpath([pwd,'/RBF'])
addpath([pwd,'/gmsh'])

%% Set FE Space and load mesh
fem          =  'P1';
[vertices, boundaries, elements] = msh_to_Mmesh( 'AcousticHornFine', 2);

%% Solve nonaffine FOM for a given configuration
% param_test = [900  0.02       0.01      0.02     0.03];
% vertices   =  RBF_DeformGeometry(param_test(2:end),vertices);
% [U, FE_SPACE, MESH, DATA]  = Elliptic_Solver(elements, vertices, boundaries, fem, 'horn_data', param_test);
% ADR_export_solution(2, real(U(1:MESH.numVertices)), MESH.vertices, MESH.elements(1:3,:), 'TestSolution');

%% Generate Affine FOM by DEIM
[ FOM ]      = build_basicFOM(elements, vertices, boundaries, fem, 'horn_data');
FOM.u_D      = @(x,mu) [];
save FOM FOM;

%% Build POD-based ROM
mu_train_Dimension = 250; 
mu_cube            = lhsdesign(mu_train_Dimension, FOM.P); % normalized design
mu_train           = bsxfun(@plus,FOM.mu_min,bsxfun(@times,mu_cube,(FOM.mu_max-FOM.mu_min)));
%[mu_train]         = FullFactorial_ParameterSpace(FOM.P, FOM.mu_min, FOM.mu_max, FOM.mu_bar, [5 2 2 2 2]);
%mu_train_Dimension = size(mu_train,1);

%% collect solution and system snapshots
ndf       = FOM.FE_SPACE.numDof; % number of DOFs
S_Au      = zeros(ndf,mu_train_Dimension);
S_f       = zeros(ndf,mu_train_Dimension);
S_u       = zeros(ndf,mu_train_Dimension);

MESH_tmp            =  FOM.MESH;
MESH_tmp.vertices   =  RBF_DeformGeometry(FOM.mu_bar(FOM.DATA.shape_param),FOM.MESH.vertices);
[MESH_tmp.jac, MESH_tmp.invjac, MESH_tmp.h] = geotrasf(MESH_tmp.dim, MESH_tmp.vertices, MESH_tmp.elements);
DATA_tmp            =  FOM.DATA;
DATA_tmp.param      =  FOM.mu_bar;
[A, F]    =  ADR_Assembler(MESH_tmp, DATA_tmp, FOM.FE_SPACE);
[A_0]     =  ApplyBC(A, F, FOM.FE_SPACE, MESH_tmp, DATA_tmp);

for i = 1 : mu_train_Dimension
   
      MESH_tmp            =  FOM.MESH;
      MESH_tmp.vertices   =  RBF_DeformGeometry(mu_train(i,FOM.DATA.shape_param),FOM.MESH.vertices);
      [MESH_tmp.jac, MESH_tmp.invjac, MESH_tmp.h] = geotrasf(MESH_tmp.dim, MESH_tmp.vertices, MESH_tmp.elements);
      DATA_tmp            =  FOM.DATA;
      DATA_tmp.param      =  mu_train(i,:);
      
      [A, F]              =  ADR_Assembler(MESH_tmp, DATA_tmp, FOM.FE_SPACE);
      [A_in, F_in]        =  ADR_ApplyBC(A, F, FOM.FE_SPACE, MESH_tmp, DATA_tmp);
      u = A_in  \ F_in;
      
      S_Au(:,i)          = (A_in - A_0) * u;
      S_f(:,i)           = F_in;
      S_u(:,i)           = u;
      
      fprintf('\nAssembled Snapshot Matrix and Vector %d of %d', i, mu_train_Dimension);
end

%% Reduction
fig_folder = 'Figures/DEIM/';
[~,~,~] =  mkdir(fig_folder);

tolPOD = [1e-4 1e-4 1e-4];

% V is the RB trial and test space
[V,    ~,    Sigma_sol]     = VPOD_basis_computation(S_u, FOM.Xnorm, tolPOD(3), 0);
figure
semilogy(Sigma_sol,'-or');
grid on
title('Solution snapshots eigenvalues')

[U_Au, ~, Sigma_matrix]  = VPOD_basis_computation(S_Au, FOM.Xnorm, tolPOD(1), 0);%svds(S_matrix,length(mu_train));
figure
semilogy(Sigma_matrix,'-ob');
grid on
title('Matrix snapshots eigenvalues')

[U_rhs,    ~,    Sigma_rhs]     = VPOD_basis_computation(S_f, FOM.Xnorm, tolPOD(2), 0);%svd(S_rhs,0);
figure
semilogy(Sigma_rhs,'-ok');
grid on
title('Rhs snapshots eigenvalues')

N = size(V,2);

% run DEIM over matrix and rhs basis
[IDEIM_m, P_m]     = DEIM(U_Au);
[IDEIM_rhs, P_rhs] = DEIM(U_rhs);

% build and store reduced matrices
Qa =  size(U_Au,2);
Qf =  size(U_rhs,2);

A_N0    = V'*(A_0*V);
A_N_VPHI = V' * U_Au;
F_N    = V'*U_rhs;
 

% find elements of the "reduced mesh" corresponding to DEIM indices
[ ~, node_to_element, node_to_boundary ] = compute_adjacency_elements(FOM.MESH.vertices, ...
    FOM.MESH.elements, FOM.MESH.dim, FOM.MESH.boundaries); 

[IDEIM_mA_elem, ~, IDEIM_mA_bound ]     = DEIM_Index_to_Elements('rhs', IDEIM_m,   ndf, node_to_element, ...
    node_to_boundary, FOM.MESH.internal_dof, FOM.MESH.numNodes);

[IDEIM_rhs_elem, ~, IDEIM_rhs_bound]    = DEIM_Index_to_Elements('rhs', IDEIM_rhs, ndf, node_to_element, ...
    node_to_boundary, FOM.MESH.internal_dof, FOM.MESH.numNodes);

IDEIM_all_elem       = unique([IDEIM_mA_elem  IDEIM_rhs_elem]);

ADR_export_solution(FOM.MESH.dim, ones(FOM.MESH.numVertices,1), FOM.MESH.vertices, FOM.MESH.elements(1:3,:), [fig_folder,'Reference_Mesh']);
ADR_export_solution(FOM.MESH.dim, ones(FOM.MESH.numVertices,1), FOM.MESH.vertices, FOM.MESH.elements(1:3,IDEIM_all_elem), [fig_folder,'Reduced_Mesh']);


%% ONLINE

% build testing set
mu_test_Dimension = 40;% [1];
mu_cube           = lhsdesign(mu_test_Dimension, FOM.P); % normalized design
mu_test           = bsxfun(@plus,FOM.mu_min,bsxfun(@times,mu_cube,(FOM.mu_max-FOM.mu_min)));%mu_train;%FOM.mu_bar;%%bsxfun(@plus,FOM.mu_min,bsxfun(@times,mu_cube,(FOM.mu_max-FOM.mu_min)));

error_sol = zeros(mu_test_Dimension,1);

I = FOM.MESH.internal_dof;
% solve the Reduced and Full Order Model on the testing set to evaluate the
% error
for i = 1 : size(mu_test,1)
    
    MESH_tmp            =  FOM.MESH;
    MESH_tmp.vertices   =  RBF_DeformGeometry(mu_test(i,FOM.DATA.shape_param),FOM.MESH.vertices);
    MESH_tmp.elements   =  FOM.MESH.elements(:,IDEIM_all_elem);
    MESH_tmp.numElem    = size(MESH_tmp.elements,2);

    [MESH_tmp.jac, MESH_tmp.invjac, MESH_tmp.h] = geotrasf(MESH_tmp.dim, MESH_tmp.vertices, MESH_tmp.elements);
    DATA_tmp            =  FOM.DATA;
    DATA_tmp.param      =  mu_test(i,:);
     
    [A, F]              =  ADR_Assembler(MESH_tmp, DATA_tmp, FOM.FE_SPACE);
    [A_in, F_in]        =  ADR_ApplyBC(A, F, FOM.FE_SPACE, MESH_tmp, DATA_tmp);
     
    A_d = A_in(IDEIM_m,:)  - A_0(IDEIM_m,:);
    A_r = A_N_VPHI * ( U_Au(IDEIM_m,:) \ (A_d * V) ) + A_N0;
     
    alpha_rhs         = U_rhs(IDEIM_rhs,:)\F_in(IDEIM_rhs);
    F_r = F_N * alpha_rhs;
    
    % solve reduced model
    U_r = A_r\F_r;
     
    uNh                            = zeros(FOM.MESH.numNodes,1);
    uNh(FOM.MESH.internal_dof)     = V*U_r;
    uNh(FOM.MESH.Dirichlet_dof)    = FOM.u_D(FOM.MESH.nodes(:,FOM.MESH.Dirichlet_dof), mu_test(i,:));
    ADR_export_solution(2, real(uNh(1:FOM.MESH.numVertices)), MESH_tmp.vertices, FOM.MESH.elements(1:3,:), 'Figures/Horn_RBsol');

    % solve Full Order Model
    [Uh]  = Elliptic_Solver(2, elements, MESH_tmp.vertices, boundaries, fem, 'horn_data', mu_test(i,:));
    
    % compute errors
    error_sol(i)      = sqrt(abs((Uh(I) - V*U_r)'*FOM.Xnorm*(Uh(I) - V*U_r))) / sqrt(abs(Uh(I)'*FOM.Xnorm*Uh(I)));
       
end
mean(error_sol)
