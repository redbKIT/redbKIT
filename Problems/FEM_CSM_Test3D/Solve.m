%clc
clear all

dim      =  3;
fem      =  'P2';

%% load P1 mesh
[vertices, boundaries, elements] = msh_to_Mmesh('beam3D_Coarse', dim);

%% Solve
[U, FE_SPACE, MESH, DATA] = CSM_Solver(dim, elements, vertices, boundaries, fem, 'datafile');