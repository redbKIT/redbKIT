clear all
clc
close all

%% Set FE Space and load mesh
fem          =  'P1';
[vertices, boundaries, elements] = msh_to_Mmesh( '../mesh/obstacle_Coarse', 2);


%% Generate Affine FOM by (M)DEIM
tolPOD_MDEIM            =  [1e-6 1e-6];
[ FOM_SEMMT ]           =  build_affineFOM(elements, vertices, boundaries, fem, 'SEMMT_data', tolPOD_MDEIM);
FOM_SEMMT.u_D{1}        =  @(x,mu)(FOM_SEMMT.DATA.bcDir{1}(x(1,:),x(2,:),[],mu));
FOM_SEMMT.u_D{2}        =  @(x,mu)(FOM_SEMMT.DATA.bcDir{2}(x(1,:),x(2,:),[],mu));

save FOM_SEMMT FOM_SEMMT;

%% Build POD-based ROM
mu_train_Dimension = 50; 
mu_cube            = lhsdesign(mu_train_Dimension, FOM_SEMMT.P); % normalized design
mu_train           = bsxfun(@plus,FOM_SEMMT.mu_min,bsxfun(@times,mu_cube,(FOM_SEMMT.mu_max-FOM_SEMMT.mu_min)));
%[mu_train]        = FullFactorial_ParameterSpace(FOM.P, FOM.mu_min, FOM.mu_max, FOM.mu_ref, [12 2 2 2 2]);

tolPOD             = 1e-6;
method             = 'Galerkin';%'LeastSquares'
[ ROM_SEMMT ]      = build_PODbased_ROM(FOM_SEMMT, mu_train, tolPOD, method);

ROM_SEMMT.u_D{1} = @(x,mu)(deform_boundary(x(1,:),x(2,:),[],[1 mu]));
ROM_SEMMT.u_D{2} = @(x,mu)(deform_boundary(x(1,:),x(2,:),[],[2 mu]));
save ROM_SEMMT ROM_SEMMT;

%% Solve ROM online, evaluate error estimate and visualize solution
mu         =  FOM_SEMMT.mu_max;
t_RB = tic;
[uN, uNh]  = solve_RBsystem(ROM_SEMMT, mu);
t_RB = toc(t_RB);
fprintf('\n -- RB problem solved in %2.2e s -- \n', t_RB)
%[ deltaN ] = error_estimate(ROM_SEMMT, uN, mu);

CSM_export_solution(2, uNh, ROM_SEMMT.MESH.vertices+[uNh(1:ROM_SEMMT.MESH.numVertices)';uNh(1+ROM_SEMMT.MESH.numNodes:ROM_SEMMT.MESH.numNodes+ROM_SEMMT.MESH.numVertices)'],...
    ROM_SEMMT.MESH.elements, ROM_SEMMT.MESH.numNodes, 'SEMMT_RB_Deformed');

CSM_export_solution(2, uNh, ROM_SEMMT.MESH.vertices,...
    ROM_SEMMT.MESH.elements, ROM_SEMMT.MESH.numNodes, 'SEMMT_RB_reference');

%% Perform Error analysis
%[ROM_SEMMT.Cqq, ROM_SEMMT.dqq, ROM_SEMMT.Eqq] = offline_residual(FOM_SEMMT, ROM_SEMMT.V);
Ntest_sample  = 15;
%[ Error ]     = error_analysis(FOM_SEMMT, ROM_SEMMT, Ntest_sample);
[ Error ]     = CSM_error_analysis_MDEIM(FOM_SEMMT, ROM_SEMMT, Ntest_sample, 'SEMMT_data', elements, vertices, boundaries);