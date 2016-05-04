clc
clear all

dim      =  2;

%% load P1 mesh
[vertices, boundaries, elements] = msh_to_Mmesh('channelexp2', dim);


%% P2-P1 approximation
fem        = {'P2', 'P1'};

[U, MESH, DATA] = NS_Solver(dim, elements, vertices, boundaries, fem, 'datafile', [], 'SolP2');