function [] = test()
%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

%% Set FE Space and load mesh 
fem      =  'P1';
[vertices, boundaries, elements] = msh_to_Mmesh( 'Device_mesh', 2);


%% Generate Affine Full-Order Model
[ FOM ] = build_affineFOM( elements, vertices, boundaries, fem, 'Device_Heat_data' );

%% Build Approximation of the Stability Factor (by RBF interpolation)
FOM.stabFactor.mu_interp_index   = [1 2 3];
FOM.stabFactor.interp_step       = [2 2 2];
FOM.stabFactor.rbf_parfor        = 1;
FOM.stabFactor.inf_sup           = 1;
[FOM]                            = RBF_OfflineInterpolation(FOM);

%% Generate Reduced-Order Model by POD

% Define sample grid where to compute POD snapshots (here LHS design)
mu_train_Dimension = 10; 
mu_cube            = lhsdesign(mu_train_Dimension, FOM.P); % normalized design
mu_train           = bsxfun(@plus,FOM.mu_min,bsxfun(@times,mu_cube,(FOM.mu_max-FOM.mu_min)));

tolPOD             = 1e-2;
method             = 'Galerkin';%'LeastSquares'
[ ROM ]            = build_PODbased_ROM(FOM, mu_train, tolPOD, method);

%% Generate Reduced-Order Model by Greedy Algorithm

% % Define sample grid where to evaluate error estimate (here LHS design)
mu_train_Dimension = 50;  
mu_cube            = lhsdesign(mu_train_Dimension, FOM.P); % normalized design
mu_train           = bsxfun(@plus,FOM.mu_min,bsxfun(@times,mu_cube,(FOM.mu_max-FOM.mu_min)));

tolGREEDY          = 1e-2;
method             = 'Galerkin';%'LeastSquares'
Nmax               = 5;
mu_1               = [0.0 20 20];
[ ROM ]            = build_GREEDYbased_ROM(FOM, mu_train, tolGREEDY, mu_1, Nmax, method);

%% Solve ROM online, evaluate error estimate and visualize solution
mu         =  [0.15   12   34 ];
t_RB = tic;
[uN, uNh]  = solve_RBsystem(ROM, mu);
t_RB = toc(t_RB);
fprintf('\n -- RB problem solve in %2.2e s -- \n', t_RB)
[ deltaN ] = error_estimate(ROM, uN, mu);

% get deformed mesh
def_vertices  = Device_Heat_DefGeo(ROM.MESH.vertices, mu);

%% Perform Error analysis
Ntest_sample  = 10;
[ Error ]     = error_analysis(FOM, ROM, Ntest_sample);

return
