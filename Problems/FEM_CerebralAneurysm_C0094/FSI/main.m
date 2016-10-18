%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri@epfl.ch>

clear all
clc

[~,~,~] = mkdir('Figures');
[~,~,~] = mkdir('Results');

dim      =  3;

%% Load F mesh
[vertices, boundaries, elements] = msh_to_Mmesh('../mesh/FluidVeryCoarse', dim);

%% Solve Fluid
[U0_Fluid, MESH, DATA] = NSt_Solver(dim, elements, vertices, boundaries, {'P1','P1'}, 'datafile_CFD_Steady', [], 'Figures/AneurysmC0094_PreStressFluid_');

save U0_Fluid U0_Fluid;

%% Solve Solid
[mshS.vertices, mshS.boundaries, mshS.elements, mshS.rings] = msh_to_Mmesh('../mesh/SolidVeryCoarse', dim);
[mshF.vertices, mshF.boundaries, mshF.elements, mshF.rings] = msh_to_Mmesh('../mesh/FluidVeryCoarse', dim);

[R_P, J_P] = FSI_PrestressSolver(dim, mshF, mshS, {'P1','P1'}, 'P1', 'datafile_CFD', 'datafile_CSM_Prestress', [], 'Figures/AneurysmC0094_PreStressSolid_');

save R_P R_P;
save J_P J_P;

%% Solve FSI
FSIt_Solver(dim, mshF, mshS, {'P1','P1'}, 'P1', 'datafile_CFD', 'datafile_CSM', [], 'Figures/AneurysmC0094_FSIPrestress_');
