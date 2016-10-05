%   For reference, see Section 3.8, 6.6 and 7.2 of
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

dim = 3;

%% Solve Steady Navier-Stokes equations using SUPG stabilized P1-P1 Finite Elements to compute the advection field
femNS      = {'P1', 'P1'};
[~,~,~]    = mkdir('Figures');
[vertices, boundaries, elements] = msh_to_Mmesh( 'mesh/mixer_gmsh_3D_sym', dim);

Advection_Field = NS_Solver(dim, elements, vertices, boundaries, femNS, 'mixerNS_data', [], 'Figures/AdvectionField');

%% Set FE Space and load mesh 
fem      =  'P1';

%% Generate Affine Full-Order Model
[ FOM ] = build_affineFOM( elements, vertices, boundaries, fem, 'mixerADR_data', Advection_Field );

%% Build Approximation of the Stability Factor (by RBF interpolation)
FOM.stabFactor.mu_interp_index   = [4];
FOM.stabFactor.interp_step       = [24];
FOM.stabFactor.rbf_parfor        = 1;
FOM.stabFactor.inf_sup           = 0;
[FOM]                            = RBF_OfflineInterpolation(FOM);

%% Generate Reduced-Order Model by POD

% Define sample grid where to compute POD snapshots (here LHS design)
mu_train_Dimension = 100; 
mu_cube            = lhsdesign(mu_train_Dimension, FOM.P); % normalized design
mu_train           = bsxfun(@plus,FOM.mu_min,bsxfun(@times,mu_cube,(FOM.mu_max-FOM.mu_min)));

tolPOD             = 1e-5;
method             = 'Galerkin';%'LeastSquares'
[ ROM ]            = build_PODbased_ROM(FOM, mu_train, tolPOD, method);

%% Generate Reduced-Order Model by Greedy Algorithm using n_train = 100

% Define sample grid where to evaluate error estimate (here LHS design)
mu_train_Dimension = 100;  
mu_cube            = lhsdesign(mu_train_Dimension, FOM.P); % normalized design
mu_train           = bsxfun(@plus,FOM.mu_min,bsxfun(@times,mu_cube,(FOM.mu_max-FOM.mu_min)));

tolGREEDY          = 1e-4;
method             = 'Galerkin';%'LeastSquares'
Nmax               = 50;
mu_1               = [6 6 6 300];
t_g100             = tic;
[ ROM_g100 ]       = build_GREEDYbased_ROM(FOM, mu_train, tolGREEDY, mu_1, Nmax, method);
ROM_g100.offline_time = toc(t_g100);
save ROM_g100 ROM_g100;


%% Generate Reduced-Order Model by Greedy Algorithm using n_train = 1000

% Define sample grid where to evaluate error estimate (here LHS design)
mu_train_Dimension = 1000;  
mu_cube            = lhsdesign(mu_train_Dimension, FOM.P); % normalized design
mu_train           = bsxfun(@plus,FOM.mu_min,bsxfun(@times,mu_cube,(FOM.mu_max-FOM.mu_min)));

tolGREEDY          = 1e-4;
method             = 'Galerkin';%'LeastSquares'
Nmax               = 50;
mu_1               = [6 6 6 300];
t_g1000             = tic;
[ ROM_g1000 ]      = build_GREEDYbased_ROM(FOM, mu_train, tolGREEDY, mu_1, Nmax, method);
ROM_g1000.offline_time = toc(t_g1000);
save ROM_g1000 ROM_g1000;

%% Generate Reduced-Order Model by Greedy Algorithm using n_train = 10000

% Define sample grid where to evaluate error estimate (here LHS design)
mu_train_Dimension = 10000;  
mu_cube            = lhsdesign(mu_train_Dimension, FOM.P); % normalized design
mu_train           = bsxfun(@plus,FOM.mu_min,bsxfun(@times,mu_cube,(FOM.mu_max-FOM.mu_min)));

tolGREEDY          = 1e-4;
method             = 'Galerkin';%'LeastSquares'
Nmax               = 50;
mu_1               = [6 6 6 300];
t_g10000           = tic;
[ ROM_g10000 ]     = build_GREEDYbased_ROM(FOM, mu_train, tolGREEDY, mu_1, Nmax, method);
ROM_g10000.offline_time = toc(t_g10000);
save ROM_g10000 ROM_g10000;

%% Print Greedy Times
ROM_g100.offline_time
ROM_g1000.offline_time
ROM_g10000.offline_time

%% Compare Greedy Convergence Histories
handle = figure
semilogy(ROM_g100.delta_Max, '-ks', 'linewidth',2)
hold on
semilogy(ROM_g1000.delta_Max, '--k', 'linewidth',2)
semilogy(ROM_g10000.delta_Max, '-ko', 'linewidth',2)
grid on
legend('n_{train} = 10^2', 'n_{train} = 10^3', 'n_{train} = 10^4')
xlabel('N')
ylabel('\max_{\mu \in \Xi_{\text{train}}} \Delta_N(\mu)')
set(findall(handle,'-property','FontSize'),'FontSize',14)

%% Solve ROM online, evaluate error estimate and visualize solution
mu         =  [8 2 4 550];
t_RB = tic;
[uN, uNh]  = solve_RBsystem(ROM_g10000, mu);
t_RB = toc(t_RB);
fprintf('\n -- RB problem solved in %2.2e s -- \n', t_RB)
[ deltaN ] = error_estimate(ROM_g10000, uN, mu);
% 
% [~,~,~] = mkdir('Figures');
% ADR_export_solution(3, uNh, ROM.MESH.vertices, ROM.MESH.elements, 'Figures/SOL_RB');

%% Perform Error analysis
Ntest_sample  = 350;
[ Error ]     = error_analysis(FOM, ROM_g1000, Ntest_sample);

% To check exponential convergence
val     = Error.X_error_rel_av;
N       = [1 : length(val)]';
alpha   = log(val(1:end-1)./val(2:end))./(N(2:end)-N(1:end-1));
alpha_m = mean(alpha)
C       = val(1) / exp(-alpha_m);
hold on
semilogy(N, C*exp(-alpha_m*N), '--g', 'LineWidth', 3);
