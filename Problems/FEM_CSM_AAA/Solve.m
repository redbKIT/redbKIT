clc
clear all

dim      =  3;
fem      =  'P1';

%% load P1 mesh
[vertices, boundaries, elements] = msh_to_Mmesh('AAA_mesh', dim);

%% Solve
[~,~,~] = mkdir('Figures');
[U, FE_SPACE, MESH, DATA] = CSM_Solver(dim, elements, vertices, boundaries, fem, 'datafile', [], 'Figures/SolFEMP1_AAA');