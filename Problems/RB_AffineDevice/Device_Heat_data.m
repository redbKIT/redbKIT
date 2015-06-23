% Source term
data.force = @(x, y, t, param)( 1 + 0.*x.*y);

% Dirichlet
data.bcDir = @(x, y, t, param)( 0 + 0.*x.*y);    

% Neumann
data.bcNeu = @(x, y, t, param)(0.*x.*y);

% Robin
data.bcRob_alpha    = @(x,y,t,param)(0.*x);
data.bcRob_gamma    = @(x,y,t,param)(0.*x);
data.bcRob_fun      = @(x,y,t,param)(0.*x);

                                    
% BC flag
data.flag_dirichlet = [1];
data.flag_neumann   = [2 3 4];
data.flag_robin     = [];

% diffusion
data.diffusion = @(x, y, t, param)( 1 + 0.*x.*y);

% transport term
data.transport{1} = @(x, y, t, param)(  0.*x );
data.transport{2} = @(x, y, t, param)(  (x-2/3).*(1-x) + 0.*x.*y);

% reaction
data.reaction = @(x, y, t, param)(0 + 0.*x.*y);


