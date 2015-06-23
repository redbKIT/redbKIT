% Source term
data.force = @(x, y, t, param)(0.*x.*y);%(  1*exp(-((x-param(2)).^2 + (y-param(3)).^2)/(2*0.05^2) ) + 0.*x.*y );

% Dirichlet
data.bcDir = @(x, y, t, param)( 0 + 0.*x.*y);    

% Neumann
data.bcNeu = @(x, y, t, param)(0.*x.*y);

% Robin
data.bcRob_alpha    = @(x,y,t,param)(0.*x);
data.bcRob_gamma    = @(x,y,t,param)(0.*x);
data.bcRob_fun      = @(x,y,t,param)(0.*x);

                                    
% BC flag
data.flag_dirichlet = [];
data.flag_neumann   = [1 2 3 4];
data.flag_robin     = [];

% diffusion
data.diffusion    = @(x, y, t, param)( 3*1e-2 + 0.*x.*y);

% transport term
data.transport{1} = @(x, y, t, param)(  1  + 0.*x.*y);
data.transport{2} = @(x, y, t, param)(  1  + 0.*x.*y);

% reaction
data.reaction = @(x, y, t, param)(0.5 + 0.*x.*y);


