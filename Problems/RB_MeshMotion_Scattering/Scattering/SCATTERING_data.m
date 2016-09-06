%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

% Source term
data.force = @(x, y, t, param)(0.*x.*y);

% Dirichlet
data.bcDir = @(x, y, t, param)( -exp( -complex(0,1)*param(4)* ( cos(param(5))*x + sin(param(5)) *y ) ) + 0.*x.*y);    

% Neumann
data.bcNeu = @(x, y, t, param)(0.*x.*y);

data.bcRob_fun = @(x, y, t, param)(0.*x.*y);
                                    
data.bcRob_alpha = @(x, y, t, param)( -complex(0,1)*param(4) + 0.*x.*y );
                                                                                                     
% BC flag
data.flag_dirichlet = [2];
data.flag_neumann   = [];
data.flag_robin     = [1];

% diffusion
data.diffusion = @(x, y, t, param)(1 + 0.*x.*y);

% transport term
data.transport{1} = @(x, y, t, param)(0.*x.*y);
data.transport{2} = @(x, y, t, param)(0.*x.*y);

% reaction
data.reaction = @(x, y, t, param)(-param(4)^2 + 0.*x.*y);
