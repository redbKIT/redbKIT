%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

clc
clear all

fem      =  'P1';
dim      = 2;

%% build P1 mesh
[vertices, boundaries, elements] = msh_to_Mmesh('hole', dim);

%% Solve using Solid-Extension Mesh Moving Technique (SEMMT)
d_x              = 0.1;
d_y              = 0.2;
Stiffening_Power = 0.6;

[U, FE_SPACE, MESH, DATA] = CSM_Solver(dim, elements, vertices, boundaries, fem, 'SEMMT_data', [d_x d_y Stiffening_Power],   'SEMMT_Solution');
CSM_export_solution(dim, U, MESH.vertices+[U(1:MESH.numVertices)';U(1+MESH.numNodes:end)'], MESH.elements, MESH.numVertices, 'SEMMT_Deformed');



%% Solve using Harmonic-Extension (HE) Mesh Moving Technique
 
[dX, FE_SPACE, MESH, DATA]  = Elliptic_Solver(dim, elements, vertices, boundaries, fem, 'HE_data', [1 d_x d_y]);
[dY, FE_SPACE, MESH, DATA]  = Elliptic_Solver(dim, elements, vertices, boundaries, fem, 'HE_data', [2 d_x d_y]);

CSM_export_solution(dim, [dX;dY], [dX';dY'], MESH.elements, MESH.numVertices, 'HT_Solution');
CSM_export_solution(dim, [dX;dY], MESH.vertices+[dX';dY'], MESH.elements, MESH.numVertices, 'HT_Deformed');
