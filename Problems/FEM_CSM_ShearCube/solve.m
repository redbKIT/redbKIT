%FEM_CSM_ShearCube implements the test proposed in the following preprint (Sect. 6.1):
%
%   D. Bonomi, A. Manzoni, A. Quarteroni. A matrix discrete empirical interpolation
%   method for the efficient model reduction of parametrized nonlinear PDEs: application to non-
%   linear elasticity problems, 2016. Available as MATHICSE report Nr.
%   21.2016, see http://mathicse.epfl.ch/page-68906-en.html.

clc
clear all

dim      =  3;
fem      =  'P2';

%% load P1 mesh
[vertices, boundaries, elements] = msh_to_Mmesh('Cube', dim);

%% Solve
[U, FE_SPACE, MESH, DATA] = CSM_Solver(dim, elements, vertices, boundaries, fem, 'datafile', [], 'Displacement_ShearCube');