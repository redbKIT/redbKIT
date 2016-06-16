%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

clear all
clc

%% Set FE Space and load mesh
fem          =  'P1';
[vertices, boundaries, elements] = msh_to_Mmesh( '../mesh/AcousticHorn_Coarse', 2);

%% Solve nonaffine FOM for a given configuration
% param_test = [900  0.02       0.01      0.02     0.03];
% def_vertices   =  RBF_MeshDeformation(param_test(2:end),vertices);
% [U, FE_SPACE, MESH, DATA]  = Elliptic_Solver(2, elements, def_vertices, boundaries, fem, 'horn_data', param_test);
% ADR_export_solution(2, real(U(1:MESH.numVertices)), def_vertices, MESH.elements(1:3,:), 'TestSolution');

%% Generate Affine FOM by (M)DEIM
tolPOD_MDEIM = [4 95];% tolerances for POD on RHS and Matrix snapshots, respectively
[ FOM ]      = build_affineFOM(elements, vertices, boundaries, fem, 'horn_data', tolPOD_MDEIM);
FOM.u_D      = @(x,mu) [];
save FOM FOM;

%% Build POD-based ROM
mu_train_Dimension = 200; 
mu_cube            = lhsdesign(mu_train_Dimension, FOM.P); % normalized design
mu_train           = bsxfun(@plus,FOM.mu_min,bsxfun(@times,mu_cube,(FOM.mu_max-FOM.mu_min)));

tolPOD             = 80;
method             = 'Galerkin';%'LeastSquares'
residual_computation = false;
[ ROM ]            = build_PODbased_ROM(FOM, mu_train, tolPOD, method, [] , residual_computation);

%% Solve ROM online, evaluate error estimate and visualize solution
mu         =  [900  0.0       0.0      0.0     0.0];
t_RB = tic;
[uN, uNh]  = solve_RBsystem(ROM, mu);
t_RB = toc(t_RB);
fprintf('\n -- RB problem solve in %2.2e s -- \n', t_RB)

% get deformed mesh
def_vertices  = RBF_MeshDeformation(mu(ROM.DATA.shape_param), ROM.MESH.vertices);

% visualize
figure
pdeplot(def_vertices,[],ROM.MESH.elements(1:3,:),'xydata',real(uNh(1:ROM.MESH.numVertices)),'xystyle','interp',...
    'zdata',real(uNh(1:ROM.MESH.numVertices)),'zstyle','continuous',...
    'colorbar','on','mesh','off');
colormap(jet);
lighting phong

% export solution to vtk for visualization in paraview
[~,~,~] = mkdir('Figures');
ADR_export_solution(2, real(uNh(1:ROM.MESH.numVertices)), def_vertices, ROM.MESH.elements(1:3,:), 'Figures/Horn_RBsol');