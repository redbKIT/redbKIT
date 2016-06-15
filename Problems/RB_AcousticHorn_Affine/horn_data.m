% Source term
data.force = @(x, y, t, param)(0.*x.*y);

% Dirichlet
data.bcDir = @(x, y, t, param)(0.*x.*y);    

% Neumann
data.bcNeu = @(x, y, t, param)(0.*x.*y);

% Robin
data.bcRob_fun = @(x, y, t, param)(  0.*x.*y);
                                    
data.bcRob_alpha = @(x, y, t, param)( 0.*x.*y );
                                                                                                     
% BC flag
data.flag_dirichlet = [];
data.flag_neumann   = [1 4];
data.flag_robin     = [2 3];

% diffusion
data.diffusion = @(x, y, t, param)(1 + 0.*x.*y);

% transport term
data.transport{1} = @(x, y, t, param)(0.*x.*y);
data.transport{2} = @(x, y, t, param)(0.*x.*y);

% reaction
% data.reaction = @(x, y, t, param)(-(FACTOR*param(1))^2 + 0.*x.*y);
data.reaction = @(x, y, t, param)( 1 + 0.*x.*y);
