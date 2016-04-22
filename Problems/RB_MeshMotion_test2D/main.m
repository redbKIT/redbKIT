clear all
clc
close all

%% Set FE Space and load mesh
fem          =  'P1';
[vertices, boundaries, elements] = msh_to_Mmesh( 'holeFine', 2);


%% Generate Affine FOM by (M)DEIM
tolPOD_MDEIM      = [1e-5 1e-5];
[ FOM ]           = build_affineFOM(elements, vertices, boundaries, fem, 'SEMMT_data', tolPOD_MDEIM);
FOM.u_D{1}        =  @(x,mu)(FOM.DATA.bcDir{1}(x(1,:),x(2,:),[],mu));
FOM.u_D{2}        =  @(x,mu)(FOM.DATA.bcDir{2}(x(1,:),x(2,:),[],mu));

save FOM FOM;

%% Build POD-based ROM
mu_train_Dimension = 100; 
mu_cube            = lhsdesign(mu_train_Dimension, FOM.P); % normalized design
mu_train           = bsxfun(@plus,FOM.mu_min,bsxfun(@times,mu_cube,(FOM.mu_max-FOM.mu_min)));
%[mu_train]        = FullFactorial_ParameterSpace(FOM.P, FOM.mu_min, FOM.mu_max, FOM.mu_ref, [12 2 2 2 2]);

tolPOD             = 1e-5;
method             = 'Galerkin';%'LeastSquares'
[ ROM ]            = build_PODbased_ROM(FOM, mu_train, tolPOD, method);

%% Solve ROM online, evaluate error estimate and visualize solution
mu         =  [0.06 0.0 0.0];
t_RB = tic;
[uN, uNh]  = solve_RBsystem(ROM, mu);
t_RB = toc(t_RB);
fprintf('\n -- RB problem solved in %2.2e s -- \n', t_RB)
[ deltaN ] = error_estimate(ROM, uN, mu);

STR_export_solution(2, uNh, ROM.MESH.vertices+[uNh(1:ROM.MESH.numVertices)';uNh(1+ROM.MESH.numNodes:end)'],...
    ROM.MESH.elements, ROM.MESH.numVertices, 'SEMMT_RB_Deformed');

STR_export_solution(2, uNh, ROM.MESH.vertices,...
    ROM.MESH.elements, ROM.MESH.numVertices, 'SEMMT_RB_reference');

%% Perform Error analysis
Ntest_sample  = 10;
[ Error ]     = error_analysis(FOM, ROM, Ntest_sample);