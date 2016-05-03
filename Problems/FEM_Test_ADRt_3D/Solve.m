function [] = Solve( fem )
%SOLVE shows how to solve a 3D parabolic problem with redbKIT

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

if nargin < 1 || isempty( fem )
    fem = 'P1';
end

dim      = 3;

[~,~,~] = mkdir('Figures');

%% laod P1 mesh
[vertices, boundaries, elements] = msh_to_Mmesh('SliceCube', dim);

%% Solve
[U, FE_SPACE, MESH, DATA]  = ADRt_Solver(dim, elements, vertices, boundaries, fem, 'datafile', [], ['Figures/Sol', fem,'_']);

end
