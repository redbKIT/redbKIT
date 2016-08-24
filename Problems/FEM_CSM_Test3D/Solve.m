%clc
clear all

dim      =  3;
fem      =  'P2';

%% load P1 mesh
[vertices, boundaries, elements] = msh_to_Mmesh('beam3D', dim);

%% Solve
[U, FE_SPACE, MESH, DATA] = CSM_Solver(dim, elements, vertices, boundaries, fem, 'datafile', [], 'Sol');