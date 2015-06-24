%   For reference, see Section 7.5 of
%
%   Quarteroni, Manzoni, Negri - REDUCED BASIS METHODS FOR PARTIAL
%   DIFFERENTIAL EQUATIONS. AN INTRODUCTION. Springer, 2015

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

% Source term
data.force = @(x, y, t, param)(1 + 0.*x.*y);

% Dirichlet
data.bcDir = @(x, y, t, param)(0 + 0.*x.*y);

% Neumann
data.bcNeu = @(x, y, t, param)(1 + 0.*x.*y);

% Robin
data.bcRob_alpha    = @(x,y,t,param)(0.*x);
data.bcRob_gamma    = @(x,y,t,param)(0.*x);
data.bcRob_fun      = @(x,y,t,param)(0.*x);


% BC flag
data.flag_dirichlet = [1 3];
data.flag_neumann   = [2 4];
data.flag_robin     = [];

% diffusion
data.diffusion    = @(x, y, t, param)(  1 + 0.*x.*y);

% transport term
data.transport{1} = @(x, y, t, param)(  0  + 0.*x.*y);
data.transport{2} = @(x, y, t, param)(  0  + 0.*x.*y);

% reaction
data.reaction     = @(x, y, t, param)(0 + 0.*x.*y);
