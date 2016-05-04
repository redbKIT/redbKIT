clc
clear all

dim      =  3;

%% load P1 mesh
[vertices, boundaries, elements] = msh_to_Mmesh('bbody3DSymCoarse', dim);

%% P2-P1 approximation
fem        = {'P2', 'P1'};
[~,~,~]    = mkdir('Figures');

[U, MESH, DATA] = NS_Solver(dim, elements, vertices, boundaries, fem, 'bbody3D_data', [], 'Figures/SolP2');