clc
clear all

%%
dim      =  3;

[~,~,~] = mkdir('Figures');

% load P1 mesh - several refinement levels are available
[vertices, boundaries, elements] = msh_to_Mmesh('../mesh/mesh_Fluid', dim);
vertices = 0.1 * vertices;

% specify FEM approximation
fem        = {'P1', 'P1'};

NSt_Solver(dim, elements, vertices, boundaries, fem, 'datafile', [], 'Figures/Sol_Fluid_');

