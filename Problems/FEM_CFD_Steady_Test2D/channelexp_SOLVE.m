clc
clear all

dim      =  2;

%% load P1 mesh
[vertices, boundaries, elements] = msh_to_Mmesh('channelexp2', dim);


%% P1-P1 approximation with Dohrman-Bochev stabilization
fem        = {'P2', 'P1'};

[U, MESH, DATA] = NSsteadySolver(dim, elements, vertices, boundaries, fem, 'datafile', [], 'SolP2');