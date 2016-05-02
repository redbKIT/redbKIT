function [] = Solve( fem )
%SOLVE shows how to solve a parabolic problem with redbKIT

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

if nargin < 1 || isempty( fem )
    fem = 'P1';
end

dim      = 2;

[~,~,~] = mkdir('Figures');

%% build P1 mesh
[vertices, boundaries, elements] = initmesh('mesh_square','Jiggle','minimum','Hgrad',1.01,'Hmax',0.2);
[vertices, boundaries, elements] = refinemesh('mesh_square',vertices, boundaries, elements);
[vertices, boundaries, elements] = refinemesh('mesh_square',vertices, boundaries, elements);
[vertices, boundaries, elements] = refinemesh('mesh_square',vertices, boundaries, elements);

%% Solve
[U, FE_SPACE, MESH, DATA]  = ADRt_Solver(dim, elements, vertices, boundaries, fem, 'datafile', [], ['Figures/Sol', fem,'_']);

end
