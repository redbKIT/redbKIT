%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri@epfl.ch>

clc
clear all

%%
dim      =  3;

[~,~,~] = mkdir('Figures');
[~,~,~] = mkdir('Results');

% load mesh
%[vertices, boundaries, elements] = msh_to_Mmesh('../mesh/FluidBL', dim);
[vertices, boundaries, elements] = msh_to_Mmesh('../mesh/FluidCoarse', dim);

vertices = 0.1 * vertices; % vertices coordinates from mm to cm

% specify FEM approximation
fem        = {'P1', 'P1'};

% solve
NSt_Solver(dim, elements, vertices, boundaries, fem, 'datafile', [], 'Figures/Sol_FluidCoarse_');