%DATAFILE for linear elasticity problem

% Source term
data.force{1} = @(x, y, z, t, param)(0.*x.*y.*z);
data.force{2} = @(x, y, z, t, param)(0.*x.*y.*z);
data.force{3} = @(x, y, z, t, param)(-0.05*1 + 0.*x.*y.*z);

% Dirichlet
data.bcDir{1} = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcDir{2} = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcDir{3} = @(x, y, z, t, param)(0.*x.*y.*z);

% Neumann
data.bcNeu{1} = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcNeu{2} = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcNeu{3} = @(x, y, z, t, param)(0 + 0.*x.*y.*z);

% Normal Pressure
data.bcPrex   = @(x, y, z, t, param)(0*x.*y.*z);

% BC flag
data.flag_dirichlet{1} = [1];
data.flag_neumann{1}   = [];%2 3 4 5 6];
data.flag_pressure{1}  = [];%6

data.flag_dirichlet{2} = [1];
data.flag_neumann{2}   = [];%2 3 4 5 6];
data.flag_pressure{2}  = [];%6

data.flag_dirichlet{3} = [1];
data.flag_neumann{3}   = [];%2 3 4 5 6];
data.flag_pressure{3}  = [];%6

% material parameters 
data.Young   = 21.5;
data.Poisson = 0.29;
data.Density = 1;