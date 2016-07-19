function [ ROM ] = build_PODbased_ROM(FOM, sample_grid, tolPOD, method, D, residual_computation)
%BUILD_PODBASED_ROM generates a reduced model by means of POD
%
%   [ ROM ] = BUILD_PODBASED_ROM(FOM, SAMPLE_GRID, TOLPOD, METHOD, D)
%   
%   Description of inputs:
%
%       FOM is struct containing the full order model
%
%       SAMPLE_GRID is a matrix of dimension number of snaphots x FOM.P
%       containing the coordinates of the parameter points where to
%       compute the snapshots
%
%       TOLPOD can be either (i) a tolerance to choose the number of POD 
%       modes to be retained according to the relative information content 
%       criterion or (ii) directly the number of POD modes to be retained.
%
%       METHOD can be either 'Galerkin' or 'LeastSquares'
%
%       D (optional) is the matrix containing the quadrature weights for
%       the POD
%
%   The output is a struct ROM.

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

if nargin < 5
    D = [];
end

if nargin < 6 || isempty(residual_computation)
    residual_computation = true;
end

%% initialize ROM struct
if isfield(FOM, 'model')
    ROM.model   = FOM.model;
else
    FOM.model = 'ADR';
    ROM.model = 'ADR';
end

ROM.MESH   = FOM.MESH;
if isfield(FOM, 'DATA')
    ROM.DATA   = FOM.DATA;
end

if isfield(FOM, 'HyRED')
    ROM.HyRED  = FOM.HyRED;
end
ROM.Qa     = FOM.Qa;
ROM.Qf     = FOM.Qf;
ROM.P      = FOM.P;
ROM.method = method;
ROM.u_D    = FOM.u_D;
ROM.D      = D;

if isfield(FOM,'stabFactor')
    ROM.stabFactor   = FOM.stabFactor;
end

Ns = size(sample_grid,1);
S  = zeros(length(FOM.MESH.internal_dof), Ns);

%% Compute snapshots
fprintf('\n Compute snapshots ...')
parfor i = 1 : Ns

    fprintf('\n   Computing snapshot %d out of %d', i, Ns);
    uh      = solve_HFsystem(FOM, sample_grid(i,:));
    S(:, i) = uh(FOM.MESH.internal_dof);

end

%% Build basis by POD
fprintf('\n Build basis by POD ...')
[ROM.V, ROM.singular_values] = POD_basis_computation(S, FOM.Xnorm, tolPOD, D);
ROM.N = size(ROM.V,2);

%% Project system
fprintf('\n Project system ...')
[ROM.ANq, ROM.FNq] = project_System(FOM, ROM.V, method);

ROM.Xnorm = ROM.V'*(FOM.Xnorm*ROM.V);

%% Compute offline residual terms
if residual_computation
    fprintf('\n   Compute offline residual terms ...')
    [ROM.Cqq, ROM.dqq, ROM.Eqq] = offline_residual(FOM, ROM.V);
end

fprintf('\n *** Finished *** \n')

end
