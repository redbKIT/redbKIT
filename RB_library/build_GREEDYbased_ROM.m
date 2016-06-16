function [ ROM ] = build_GREEDYbased_ROM(FOM, sample_grid, tolGREEDY, mu_1, Nmax, method)
%BUILD_GREEDYBASED_ROM generates a reduced model by means of POD
%
%   [ ROM ] = BUILD_GREEDYBASED_ROM(FOM, SAMPLE_GRID, TOLGREEDY, MU_1, ...
%                                   NMAX, METHOD)
%   
%   Description of inputs:
%
%       FOM is struct containing the full order model
%
%       SAMPLE_GRID is a matrix of dimension n_train x FOM.P
%       containing the coordinates of the parameter points belonging to the
%       training set
%
%       TOLGREEDY tolerance on the (relative) error estimate
% 
%       MU_1 is the starting parameter point
%
%       NMAX maximum number of iterations
%
%       METHOD can be either 'Galerkin' or 'LeastSquares'
%
%   The output is a struct ROM.

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

%% initialize ROM struct
ROM.MESH   = FOM.MESH;
ROM.Qa     = FOM.Qa;
ROM.Qf     = FOM.Qf;
ROM.P      = FOM.P;
ROM.method = method;
ROM.u_D    = FOM.u_D;

if isfield(FOM, 'model')
    ROM.model   = FOM.model;
else
    ROM.model = 'ADR';
    FOM.model   = 'ADR';
end

ROM.MESH   = FOM.MESH;
if isfield(FOM, 'DATA')
    ROM.DATA   = FOM.DATA;
end

if isfield(FOM, 'HyRED')
    ROM.HyRED  = FOM.HyRED;
end

if isfield(FOM,'stabFactor')
    ROM.stabFactor   = FOM.stabFactor;
end

%% Initialize greedy variables
ROM.N              = 0;
delta_Max          = tolGREEDY + 1;
sample_grid        = [mu_1; sample_grid];
new_mu_indx        = 1;
ROM.V              = [];

%% Start Greedy Algorithm
while ROM.N < Nmax && delta_Max > tolGREEDY
   
    ROM.N   =  ROM.N + 1;
    fprintf('\n **** Greedy iteration number %d ', ROM.N);
    
    fprintf('\n   Compute snapshot ...')
    uh      =  solve_HFsystem(FOM, sample_grid(new_mu_indx,:));
    
    zeta_n  = Gram_Schmidt_orth(ROM.V, uh(FOM.MESH.internal_dof), FOM.Xnorm);
    ROM.V   = [ROM.V zeta_n];
    ROM.Greedy_samples(ROM.N,:) = sample_grid(new_mu_indx,:);
    
    fprintf('\n   Project system ...')
    [ROM.ANq, ROM.FNq] = project_System(FOM, ROM.V, method);

    fprintf('\n   Compute offline residual terms ...')
    [ROM.Cqq, ROM.dqq, ROM.Eqq] = offline_residual(FOM, ROM.V);
    
    fprintf('\n   Evaluate error estimate ... ')
    delta_N = zeros(1,size(sample_grid,1));
    parfor i = 1 : size(sample_grid,1)
        mu              =  sample_grid(i,:);
        uN              =  solve_RBsystem(ROM, mu);
        delta_N(i)      =  error_estimate(ROM, uN, mu) / norm(uN); 
    end
    
    [delta_Max, new_mu_indx] = max(delta_N);
    fprintf(' Max Delta_N = %2.3e \n', delta_Max)
    
    ROM.delta_Max(ROM.N) = delta_Max;
    
end

ROM.Xnorm = ROM.V'*(FOM.Xnorm*ROM.V);

fprintf('\n *** Finished *** \n')

end