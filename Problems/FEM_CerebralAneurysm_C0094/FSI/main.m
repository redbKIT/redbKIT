%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri@epfl.ch>

clear all
clc

[~,~,~] = mkdir('Figures');
[~,~,~] = mkdir('Results');

dim      =  3;

%% Load meshes
[mshS.vertices, mshS.boundaries, mshS.elements, mshS.rings] = msh_to_Mmesh('../mesh/SolidCoarse', dim);
[mshF.vertices, mshF.boundaries, mshF.elements, mshF.rings] = msh_to_Mmesh('../mesh/FluidCoarse', dim);

mshS.vertices = 0.1 * mshS.vertices;
mshF.vertices = 0.1 * mshF.vertices;

%% Solve
[U, MESH, DATA] = FSIt_Solver(dim, mshF, mshS, {'P1','P1'}, 'P1', 'datafile_CFD', 'datafile_CSM', [], 'Figures/C0094_Coarse');
