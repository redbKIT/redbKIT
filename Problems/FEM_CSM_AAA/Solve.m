%clc
clear all

dim      =  3;
fem      =  'B1';

%% load P1 mesh
[vertices, boundaries, elements] = msh_to_Mmesh('AAA_meshFineFlags', dim);
%[vertices, boundaries, elements] = msh_to_Mmesh('AAA_meshMediumFlags', dim);

%% Solve
[~,~,~] = mkdir('Figures');
[U, FE_SPACE, MESH, DATA] = CSM_Solver(dim, elements, vertices, boundaries, fem, 'datafile', [], 'Figures/NeoHookean_FineB1');