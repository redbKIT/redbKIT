function [A, F, M] = Assembler_2D(MESH, DATA, FE_SPACE, OPERATOR, TC_d, TC_t, subdomain, t)
%ASSEMBLER_2D assembler for 2D ADR equations with numerical
%quadratures.
%
%   [A, F] = ASSEMBLER_2D(MESH, DATA, FE_SPACE) given MESH and 
%   FE_SPACE structs, assemble the stiffness matrix A and rhs vector F
%   of the advection-diffusion-reaction differential operator using the
%   parameters specified in the DATA struct. The resulting matrix and
%   vector doesn't take into account for Boundary Conditions.
%
%   [A, F, M] = ASSEMBLER_2D(MESH, DATA, FE_SPACE, T) as before, but also
%   returns the mass matrix M. Moreover, the problem's parameters are
%   evaluated with respect to a given time T.
%
%   [A] = Assembler_2D(MESH, DATA, FE_SPACE, 'diffusion', [i j], [], k);
%   returns the matrix corresponding to the (i,j)-th second derivative
%   (diffusion operator) assembled over the subdomain \Omega_k
%   A spatial dependent coefficient has to be defined in DATA.diffusion
%
%   [A] = Assembler_2D(MESH, DATA, FE_SPACE, 'transport', [], [i], k);
%   returns the matrix corresponding to the (i)-th first derivative 
%   (advection operator) assembled over the subdomain \Omega_k
%   A spatial dependent coefficient has to be defined in DATA.transport{i}

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri@epfl.ch> 

warning('Assembler_2D is deprecated. Use ADR_Assembler instead.')

if nargin < 4 || isempty(OPERATOR)
    OPERATOR = 'all';
end

if nargin < 5 || isempty(TC_d)
    TC_d = [10 10];% diagonal components of the diffusion operator by default
end

if nargin < 6 || isempty(TC_t)
    TC_t = 10; % all components of the transport operator by default
end

if nargin < 7
    subdomain = [];
end

if nargin < 8
    t = [];
end

if ~isfield(FE_SPACE, 'dphi_ref')
    FE_SPACE.dphi_ref(:,:,1) = FE_SPACE.dcsiphi;
    FE_SPACE.dphi_ref(:,:,2) = FE_SPACE.detaphi;
end

[A, F, M] = ADR_Assembler(MESH, DATA, FE_SPACE, OPERATOR, TC_d, TC_t, subdomain, t);

return