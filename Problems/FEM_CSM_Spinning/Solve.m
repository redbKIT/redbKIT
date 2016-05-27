%SOLVE shows how to perform a CSM dynamic simulation

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>


fem            = 'P1';
dim      =  3;

%% load P1 mesh
[vertices, boundaries, elements] = msh_to_Mmesh('cone_SuperCoarse', dim);


%% Solve
[~,~,~] = mkdir('FiguresB');
%[U, FE_SPACE, MESH, DATA] = Spinning_Solver(dim, elements, vertices, boundaries, fem, 'Spinning3D_data', [], 'Figures/SolSpinning_');
[U, FE_SPACE, MESH, DATA] = Spinning_Solver_Alpha(dim, elements, vertices, boundaries, fem, 'Spinning3D_data', [], 'FiguresB/SolSpinningP1_');