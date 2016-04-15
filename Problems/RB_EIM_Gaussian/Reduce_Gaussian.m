%   For reference, see Section 10.5 of
%
%   Quarteroni, Manzoni, Negri - REDUCED BASIS METHODS FOR PARTIAL
%   DIFFERENTIAL EQUATIONS. AN INTRODUCTION. Springer, 2015

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

clear all
clc

%% Set FE Space and load mesh 
fem      =  'P1';
[vertices, boundaries, elements] = initmesh('mesh_square','Jiggle','minimum','Hgrad',1.01,'Hmax', 0.1);
for i = 1 : 1
    [vertices, boundaries, elements] = refinemesh('mesh_square',vertices, boundaries, elements);
end

%% Generate Affine Full-Order Model
[ FOM ] = build_affineFOM( elements, vertices, boundaries, fem, 'Gaussian_data' );
save FOM FOM;

%% Build Approximation of the Stability Factor (by RBF interpolation)
FOM.stabFactor.mu_interp_index   = [3];
FOM.stabFactor.interp_step       = [8];
FOM.stabFactor.rbf_parfor        = 1;
FOM.stabFactor.inf_sup           = 1;
[FOM]                            = RBF_OfflineInterpolation(FOM);

%% Generate Reduced-Order Model by POD

% Define sample grid where to compute POD snapshots (here LHS design)
mu_train_Dimension = 150; 
mu_cube            = lhsdesign(mu_train_Dimension, FOM.P); % normalized design
mu_train           = bsxfun(@plus,FOM.mu_min,bsxfun(@times,mu_cube,(FOM.mu_max-FOM.mu_min)));

tolPOD             = 85;
method             = 'Galerkin';%'LeastSquares'
[ ROM ]            = build_PODbased_ROM(FOM, mu_train, tolPOD, method);

save ROM ROM;

%% Generate Reduced-Order Model by Greedy Algorithm

% Define sample grid where to evaluate error estimate (here LHS design)
% mu_train_Dimension = 2000;  
% mu_cube            = lhsdesign(mu_train_Dimension, FOM.P); % normalized design
% mu_train           = bsxfun(@plus,FOM.mu_min,bsxfun(@times,mu_cube,(FOM.mu_max-FOM.mu_min)));
% 
% tolGREEDY          = 1e-4;
% method             = 'Galerkin';%'LeastSquares'
% Nmax               = 60;
% mu_1               = [0.5 0.25 0];
% [ ROM ]            = build_GREEDYbased_ROM(FOM, mu_train, tolGREEDY, mu_1, Nmax, method);

%% Solve ROM online, evaluate error estimate and export solutions to VTK
[~,~,~] = mkdir('Figures');

mu         =  [0.25     0.3   230 ];
t_RB = tic;
[uN, uNh]  = solve_RBsystem(ROM, mu);
t_RB = toc(t_RB);
fprintf('\n -- RB problem solve in %2.2e s -- \n', t_RB)
[ deltaN ] = error_estimate(ROM, uN, mu);

ADR_export_solution(2, uNh(1:ROM.MESH.numVertices), vertices, ROM.MESH.elements(1:3,:), 'Figures/RBsol_S1');

mu         =  [0.65     0.22   115 ];
[uN, uNh]  = solve_RBsystem(ROM, mu);
ADR_export_solution(2, uNh(1:ROM.MESH.numVertices), vertices, ROM.MESH.elements(1:3,:), 'Figures/RBsol_S2');

mu         =  [0.5     0.17   32 ];
[uN, uNh]  = solve_RBsystem(ROM, mu);
ADR_export_solution(2, uNh(1:ROM.MESH.numVertices), vertices, ROM.MESH.elements(1:3,:), 'Figures/RBsol_S3');

mu         =  [0.75     0.33   285 ];
[uN, uNh]  = solve_RBsystem(ROM, mu);
ADR_export_solution(2, uNh(1:ROM.MESH.numVertices), vertices, ROM.MESH.elements(1:3,:), 'Figures/RBsol_S4');

%% Perform Error analysis

% Ntest_sample  = 2;
% % [ Error ]     = error_analysis(FOM, ROM, Ntest_sample);
% [ Error ]     = EIM_error_analysis(FOM, ROM, Ntest_sample, 'GaussianFOM_data');

clear Error
Ntest_sample = 100;
N_vec = [1:5:85 85];
M = [5 10 20 30];
parfor k =  1 : length(M)
    Error{k} = EIM_error_analysis(FOM, ROM, Ntest_sample, 'GaussianFOM_data', N_vec, M(k));
end

handle = figure;
for i =  1 : length(M)
    semilogy(N_vec, Error{i}.X_error_rel_av,'-o','Color',[rand(1) rand(1) rand(1)],'linewidth',2)
    legend_string{i} = ['Q_f = ',num2str(M(i))];
    hold on
end
title(['Convergence N,M. Rel. error averaged on ', num2str(Ntest_sample) ,' samples'])
legend(legend_string)
grid on
saveas(handle,'Figures/convergence_N_M','fig');
