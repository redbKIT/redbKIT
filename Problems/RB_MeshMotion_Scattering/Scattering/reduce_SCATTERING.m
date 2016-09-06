%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

clear all
clc
close all

%% Set FE Space and load mesh
fem          =  'P1';
[vertices, boundaries, elements] = msh_to_Mmesh( '../mesh/obstacle_Coarse', 2);

%% Generate Affine FOM by (M)DEIM
loaded = load('../MeshMotion/ROM_SEMMT.mat');
ROM_SEMMT = loaded.ROM_SEMMT;
ROM_SEMMT.u_D{1} = @(x,mu)(deform_boundary(x(1,:),x(2,:),[],[1 mu]));
ROM_SEMMT.u_D{2} = @(x,mu)(deform_boundary(x(1,:),x(2,:),[],[2 mu]));

tolPOD_MDEIM = [1e-7 1e-7];
[ FOM ]      = build_affineFOM(elements, vertices, boundaries, fem, 'SCATTERING_data', tolPOD_MDEIM, ROM_SEMMT);
FOM.u_D      = @(x,mu) -exp( -complex(0,1)*mu(4)* ( cos(mu(5))*x(1,:) + sin(mu(5)) *x(2,:) ) );
save FOM FOM;

%% Build POD-based ROM
mu_train_Dimension = 60; 
mu_cube            = lhsdesign(mu_train_Dimension, FOM.P); % normalized design
mu_train           = bsxfun(@plus,FOM.mu_min,bsxfun(@times,mu_cube,(FOM.mu_max-FOM.mu_min)));

tolPOD             = 20;
method             = 'Galerkin';
[ ROM ]            = build_PODbased_ROM(FOM, mu_train, tolPOD, method, [], 0);
ROM.u_D            = @(x,mu) -exp( -complex(0,1)*mu(4)* ( cos(mu(5))*x(1,:) + sin(mu(5)) *x(2,:) ) );
save ROM ROM;

%% Solve ROM online, evaluate error estimate and visualize solution
mu         =  [0.2  1.1  0.9 2.8  0];
% get deformed mesh
def_nodes  = MoveMesh(ROM.MESH.nodes, mu(ROM.DATA.shape_param), ROM_SEMMT);
def_vertices = def_nodes(:,1:ROM.MESH.numVertices );

t_RB = tic;
[uN, uNh]  = solve_RBsystem(ROM, mu);
uNh(ROM.MESH.Dirichlet_dof)    = ROM.u_D(def_nodes(:,ROM.MESH.Dirichlet_dof), mu);
%[uNh]  = solve_HFsystem(FOM, mu);
t_RB = toc(t_RB);
fprintf('\n -- RB problem solve in %2.2e s -- \n', t_RB)
%[ deltaN ] = error_estimate(ROM, uN, mu);


% visualize
figure
pdeplot(def_vertices,[],ROM.MESH.elements(1:3,:),'xydata',real(uNh(1:ROM.MESH.numVertices)),'xystyle','interp',...
    'zdata',real(uNh(1:ROM.MESH.numVertices)),'zstyle','continuous',...
    'colorbar','on','mesh','off');
colormap(jet);
lighting phong

% export solution to vtk for visualization in paraview
[~,~,~] = mkdir('Figures');
ADR_export_solution(2, real(uNh(1:ROM.MESH.numVertices)), def_vertices, ROM.MESH.elements(1:3,:), 'Figures/Scattering_RBsol');
uh = Elliptic_Solver(FOM.MESH.dim,  elements, def_vertices, boundaries, 'P1', 'SCATTERING_data', mu);
ADR_export_solution(2, real(uh(1:ROM.MESH.numVertices)), def_vertices, ROM.MESH.elements(1:3,:), 'Figures/Scattering_RBsolTrue');

%% Perform Error analysis
Ntest_sample  = 20;
[ Error ]     = ADR_error_analysis_MDEIM(FOM, ROM, Ntest_sample, 'SCATTERING_data', elements, vertices, boundaries, ROM_SEMMT);