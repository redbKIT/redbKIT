function [] = test()
%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

%% Set FE Space and load mesh 
fem      =  'P1';
[vertices, boundaries, elements] = initmesh('mesh_square','Jiggle','minimum','Hgrad',1.01,'Hmax', 0.1);
for i = 1 : 2
    [vertices, boundaries, elements] = refinemesh('mesh_square',vertices, boundaries, elements);
end

%% Generate Affine Full-Order Model
[ FOM ] = build_affineFOM( elements, vertices, boundaries, fem, 'Gaussian_data' );

%% Build Approximation of the Stability Factor (by RBF interpolation)
FOM.stabFactor.mu_interp_index   = [3];
FOM.stabFactor.interp_step       = [3];
FOM.stabFactor.rbf_parfor        = 1;
FOM.stabFactor.inf_sup           = 1;
[FOM]                            = RBF_OfflineInterpolation(FOM);

%% Generate Reduced-Order Model by POD

% Define sample grid where to compute POD snapshots (here LHS design)
mu_train_Dimension = 20; 
mu_cube            = lhsdesign(mu_train_Dimension, FOM.P); % normalized design
mu_train           = bsxfun(@plus,FOM.mu_min,bsxfun(@times,mu_cube,(FOM.mu_max-FOM.mu_min)));

tolPOD             = 8;
method             = 'Galerkin';%'LeastSquares'
[ ROM ]            = build_PODbased_ROM(FOM, mu_train, tolPOD, method);

%% Solve ROM online, evaluate error estimate 
mu         =  [0.25     0.3   230 ];
t_RB = tic;
[uN, uNh]  = solve_RBsystem(ROM, mu);
t_RB = toc(t_RB);
fprintf('\n -- RB problem solve in %2.2e s -- \n', t_RB)
[ deltaN ] = error_estimate(ROM, uN, mu);

%% Perform Error analysis

Ntest_sample = 10;
N_vec = [4 8];
M = [5 10];
for k =  1 : length(M)
    Error{k} = EIM_error_analysis(FOM, ROM, Ntest_sample, 'GaussianFOM_data', N_vec, M(k));
end

return