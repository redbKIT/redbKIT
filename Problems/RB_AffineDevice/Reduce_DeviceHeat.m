%   mu1 geometric parameter
%   mu2 Peclet in the fluid domain
%   mu3 diffusivity in domain \Omega_3
%
%   For reference, see Section 8.3 and 9.1 of
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
[vertices, boundaries, elements] = msh_to_Mmesh( 'Device_mesh', 2);


%% Generate Affine Full-Order Model
[ FOM ] = build_affineFOM( elements, vertices, boundaries, fem, 'Device_Heat_data' );

%% Build Approximation of the Stability Factor (by RBF interpolation)
FOM.stabFactor.mu_interp_index   = [1 2 3];
FOM.stabFactor.interp_step       = [3 3 3];
FOM.stabFactor.rbf_parfor        = 1;
FOM.stabFactor.inf_sup           = 1;
[FOM]                            = RBF_OfflineInterpolation(FOM);

%% Generate Reduced-Order Model by POD

% Define sample grid where to compute POD snapshots (here LHS design)
mu_train_Dimension = 100; 
mu_cube            = lhsdesign(mu_train_Dimension, FOM.P); % normalized design
mu_train           = bsxfun(@plus,FOM.mu_min,bsxfun(@times,mu_cube,(FOM.mu_max-FOM.mu_min)));

tolPOD             = 1e-4;
method             = 'Galerkin';%'LeastSquares'
[ ROM ]            = build_PODbased_ROM(FOM, mu_train, tolPOD, method);

%% Generate Reduced-Order Model by Greedy Algorithm

% % Define sample grid where to evaluate error estimate (here LHS design)
mu_train_Dimension = 1000;  
mu_cube            = lhsdesign(mu_train_Dimension, FOM.P); % normalized design
mu_train           = bsxfun(@plus,FOM.mu_min,bsxfun(@times,mu_cube,(FOM.mu_max-FOM.mu_min)));

tolGREEDY          = 1e-4;
method             = 'Galerkin';%'LeastSquares'
Nmax               = 35;
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

% visualize
figure
pdeplot(def_vertices,[],ROM.MESH.elements(1:3,:),'xydata',uNh(1:ROM.MESH.numVertices),'xystyle','interp',...
    'zdata',uNh(1:ROM.MESH.numVertices),'zstyle','continuous',...
    'colorbar','on','mesh','off');
colormap(jet);
lighting phong

% export solution to vtk for visualization in paraview
[~,~,~] = mkdir('Figures');
ADR_export_solution(2, uNh(1:ROM.MESH.numVertices), def_vertices, ROM.MESH.elements(1:3,:), 'Figures/Device_RBsol');

%% Perform Error analysis
Ntest_sample  = 100;
[ Error ]     = error_analysis(FOM, ROM, Ntest_sample);
