FACTOR     = 2*pi/340;
A          = 1;
R          = 1;

% Source term
data.force = @(x, y, t, param)(0.*x.*y);

% Dirichlet
data.bcDir = @(x, y, t, param)(0.*x.*y);    

% Neumann
data.bcNeu = @(x, y, t, param)(0.*x.*y);

% Robin
data.bcRob_fun = @(x, y, t, param)(+((x+1)<1e-3 & y <=0.05)*complex(0,1)*2*FACTOR*param(1)*A + 0.*x.*y);
                                    
data.bcRob_alpha = @(x, y, t, param)(+((x+1)<1e-3 & y <=0.05)*complex(0,1)*FACTOR*param(1) + ...
                                          (~((x+1)<1e-3 & y <=0.05))*(complex(0,1)*FACTOR*param(1) + 1/(2*R)) + 0.*x.*y );
                                                                                                     
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
data.reaction = @(x, y, t, param)(-(FACTOR*param(1))^2 + 0.*x.*y);
