clc
clear all

dim      =  3;

%% load P1 mesh
[vertices, boundaries, elements] = msh_to_Mmesh('../mesh/mesh_SolidCoarse', dim);

% mesh is in mm, we convert it to cm

vertices = 0.1 * vertices;

%% Solve
[~,~,~] = mkdir('Figures');

fem      =  'P1';
pressure = 50;%linspace(120,120,10);

U = [];
for i = 1 : length( pressure )
    
    CSM_Solver(dim, elements, vertices, boundaries, fem, 'datafile', pressure(i), ['Figures/SV2_P1_',num2str(pressure(i)),'mmHg'], U);
    
end