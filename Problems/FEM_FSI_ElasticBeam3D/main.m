%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri@epfl.ch>

clear all
clc

[~,~,~] = mkdir('Figures');
[~,~,~] = mkdir('Results');

dim      =  3;

%% Load meshes
[mshS.vertices, mshS.boundaries, mshS.elements, mshS.rings] = msh_to_Mmesh('mesh/mesh_Solid', dim);
[mshF.vertices, mshF.boundaries, mshF.elements, mshF.rings] = msh_to_Mmesh('mesh/mesh_Fluid', dim);

%% Solve
[U, MESH, DATA] = FSIt_Solver(dim, mshF, mshS, {'P1','P1'}, 'P1', 'NS_data', 'CSM_data', [], 'Figures/Benchmark_ImplP1SVK');