function [] = test(fem)

if nargin < 1 || isempty(fem)
    fem = 'P1';
end

dim      =  2;

%% build P1 mesh
[vertices, boundaries, elements] = initmesh('mesh_rectangle','Jiggle','minimum','Hgrad',1.01,'Hmax',0.04);

% refine
[vertices, boundaries, elements] = refinemesh('mesh_rectangle',vertices, boundaries, elements);

%% Solve
CSM_Solver(dim, elements, vertices, boundaries, fem, 'LinearElasticity2D_data');

end