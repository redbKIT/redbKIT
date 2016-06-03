function U = Solve(material_model, tf, fem)
%SOLVE shows how to perform a CSM dynamic simulation

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

if nargin < 1 
   
    material_model = 'Linear';
    tf             = 10;
    fem            = 'P1';
end


dim      =  2;

%% build P1 mesh
[vertices, boundaries, elements] = initmesh('mesh_rectangle','Jiggle','minimum','Hgrad',1.01,'Hmax',0.04);

% refine
[vertices, boundaries, elements] = refinemesh('mesh_rectangle',vertices, boundaries, elements);
[vertices, boundaries, elements] = refinemesh('mesh_rectangle',vertices, boundaries, elements);

%% Solve
[~,~,~] = mkdir('Figures');

param{1} = material_model;
param{2} = tf;

[U, FE_SPACE, MESH, DATA] = CSMt_Solver(dim, elements, vertices, boundaries, fem, 'datafile', param, ['Figures/Sol_',material_model,fem,'_']);


end