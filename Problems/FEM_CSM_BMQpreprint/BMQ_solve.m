clc
clear all

dim      =  3;
fem      =  'P1';

%% load P1 mesh
[vertices, boundaries, elements] = msh_to_Mmesh('cube', dim);

%% Solve
[U, FE_SPACE, MESH, DATA] = CSM_Solver(dim, elements, vertices, boundaries, fem, 'BMQ_data');