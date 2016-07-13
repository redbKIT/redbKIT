clear all
clc

[~,~,~] = mkdir('Figures');
[~,~,~] = mkdir('Results');

dim      =  2;

%% Load meshes
[mshS.vertices, mshS.boundaries, mshS.elements, mshS.rings] = msh_to_Mmesh('mesh/mesh_Solid', dim);
[mshF.vertices, mshF.boundaries, mshF.elements, mshF.rings] = msh_to_Mmesh('mesh/mesh_Fluid', dim);

%% Solve
[U, MESH, DATA] = FSIt_Solver(dim, mshF, mshS, {'P1','P1'}, 'P1', 'NS_data', 'CSM_data', [], 'Figures/Benchmark_Impl_');