clc
clear all

dim      =  3;

%% load P1 mesh
[vertices, boundaries, elements] = msh_to_Mmesh('AAA_meshMediumFlags', dim);
%[vertices, boundaries, elements] = msh_to_Mmesh('AAA_meshMediumFlags', dim);

%% Solve
[~,~,~] = mkdir('Figures');

%fem      =  'P2';
%U_P2 = CSM_Solver(dim, elements, vertices, boundaries, fem, 'datafile', [], 'Figures/NeoHookean_FineP2');


fem      =  'P2';
U_P1R = CSM_Solver(dim, elements, vertices, boundaries, fem, 'datafile', {'RaghavanVorp', 1}, 'Figures/RagavanVorp_MediumP2');
