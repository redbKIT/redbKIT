clc
clear all

dim      =  2;
fem      =  'P2';

%% build P1 mesh
[vertices, boundaries, elements] = initmesh('mesh_rectangle','Jiggle','minimum','Hgrad',1.01,'Hmax',0.04);

% refine
[vertices, boundaries, elements] = refinemesh('mesh_rectangle',vertices, boundaries, elements);
[vertices, boundaries, elements] = refinemesh('mesh_rectangle',vertices, boundaries, elements);

%% Solve
[U, FE_SPACE, MESH, DATA] = CSM_Solver(dim, elements, vertices, boundaries, fem, 'LinearElasticity2D_data', [], 'Sol');