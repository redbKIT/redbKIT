%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

[~,~,~] = mkdir( 'Figures' );

fem   = 'P1';
dim   =  3;

%% load P1 mesh
[vertices, boundaries, elements] = msh_to_Mmesh('mesh/cone_SuperCoarse', dim);

%% Solve
[~,~,~] = mkdir('Figures');
CSMt_Solver(dim, elements, vertices, boundaries, fem, 'Spinning3D_data', [], 'Figures/SolSpinningCoarse_', false);