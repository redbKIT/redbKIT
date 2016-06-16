clc
clear all

dim      =  3;
fem      =  'P1';

%% load P1 mesh
[vertices, boundaries, elements] = msh_to_Mmesh('Cube', dim);

%% Solve
[U, FE_SPACE, MESH, DATA] = CSM_Solver(dim, elements, vertices, boundaries, fem, 'datafile', [], 'Displacement_ShearCube');