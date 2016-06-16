%   For reference, see Section 8.6 and 9.2 of
%
%   Quarteroni, Manzoni, Negri - REDUCED BASIS METHODS FOR PARTIAL
%   DIFFERENTIAL EQUATIONS. AN INTRODUCTION. Springer, 2015
%
%   The high-fidelity (full-order) model is already provided in the FOM.mat
%   struct

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

clear all
clc

%% Load Affine Full-Order Model
load FOM FOM;
FOM.u_D   = @(x,mu) (FOM.u_D_vec);

%% Build Approximation of the Stability Factor (by RBF interpolation)
FOM.stabFactor.mu_interp_index   = [1 2];
FOM.stabFactor.interp_step       = [5 3];
FOM.stabFactor.fine_sample_size  = 1000;
FOM.stabFactor.rbf_parfor        = 1;
FOM.stabFactor.inf_sup           = 0;
[FOM]                            = RBF_OfflineInterpolation(FOM);

%% Generate Reduced-Order Model by POD
% Define sample grid where to compute POD snapshots (here LHS design)
mu_train_Dimension = 100; 
mu_cube            = lhsdesign(mu_train_Dimension, FOM.P); % normalized design
mu_train           = bsxfun(@plus,FOM.mu_min,bsxfun(@times,mu_cube,(FOM.mu_max-FOM.mu_min)));

tolPOD             = 15;
method             = 'Galerkin';%'LeastSquares'
[ ROMp ]           = build_PODbased_ROM(FOM, mu_train, tolPOD, method);

%% Generate Reduced-Order Model by Greedy Algorithm

% Define sample grid where to evaluate error estimate (here LHS design)
mu_train_Dimension = 5000;  
mu_cube            = lhsdesign(mu_train_Dimension, FOM.P); % normalized design
mu_train           = bsxfun(@plus,FOM.mu_min,bsxfun(@times,mu_cube,(FOM.mu_max-FOM.mu_min)));

tolGREEDY          = 1e-5;
method             = 'Galerkin';%'LeastSquares'
Nmax               = 25;
mu_1               = FOM.mu_min;
[ ROMg ]            = build_GREEDYbased_ROM(FOM, mu_train, tolGREEDY, mu_1, Nmax, method);

%% Solve ROM online, evaluate error estimate and visualize solution
ROM        =  ROMp;
mu         =  FOM.mu_max;
t_RB       = tic;
[uN, uNh]  = solve_RBsystem(ROM, mu);
%[uNh]  = solve_HFsystem(FOM, mu);%solve_RBsystem(ROM, mu);
t_RB       = toc(t_RB);
fprintf('\n -- RB problem solve in %2.2e s -- \n', t_RB)
[ deltaN ] = error_estimate(ROM, uN, mu);

[~,~,~]    = mkdir('Figures');
CSM_export_solution(3, uNh, ROM.MESH.vertices, ROM.MESH.elements, ROM.MESH.numVertices, 'Figures/SOL_RB_max');

%% Perform Error analysis
Ntest_sample  = 30;
[ Error ]     = error_analysis(FOM, ROMp, Ntest_sample);