%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

clc
clear all
[~,~,~] = mkdir('Figures_FEM');

fem      =  'P1';
dim      = 2;

%% load P1 mesh
[vertices, boundaries, elements] = msh_to_Mmesh('../mesh/obstacle_Coarse', dim);


%% Parameters
mu = [-0.4 1.2 1.2 5 pi/6];


%% Deform-Mesh Solid-Extension Mesh Moving Technique (SEMMT)
px               = mu(1);
py               = mu(2);
Stiffening_Power = mu(3);

[U, FE_SPACE, MESH, DATA] = CSM_Solver(dim, elements, vertices, boundaries, 'P1', 'SEMMT_data', [px py Stiffening_Power],   'Figures_FEM/SEMMT_Solution');
CSM_export_solution(dim, U, MESH.vertices+[U(1:MESH.numVertices)';U(1+MESH.numNodes:end)'], MESH.elements, MESH.numVertices, 'Figures_FEM/SEMMT_Deformed');

def_vertices = MESH.vertices+[U(1:MESH.numVertices)';U(1+MESH.numNodes:end)'];


%% Solve Helmoltz equations on deformed geometry
[U_h, FE_SPACE, MESH, DATA]  = Elliptic_Solver(dim, elements, def_vertices, boundaries, fem, 'SCATTERING_data', mu(:));
ADR_export_solution(2, real( U_h(1:MESH.numVertices) ), MESH.vertices, MESH.elements, 'Figures_FEM/Sol_FEM');

