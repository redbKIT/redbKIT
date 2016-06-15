function [] = test(fem)

if nargin < 1 || isempty(fem)
    fem = 'P1';
end

dim      =  3;

%% load P1 mesh
[vertices, boundaries, elements] = msh_to_Mmesh('beam3D', dim);

%% Solve
CSM_Solver(dim, elements, vertices, boundaries, fem, 'datafile');

end