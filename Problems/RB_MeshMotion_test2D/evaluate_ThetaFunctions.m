function [ theta_a, theta_f ] = evaluate_ThetaFunctions( mu, FOM, Ma, Mf )
%   Author: F. Negri (federico.negri@epfl.ch) 2015
%   Copyright (C) Federico Negri, CMCS, EPFL

if nargin < 3 
    Ma = FOM.Qa;
end

if nargin < 4 
    Mf = FOM.Qf;
end

n_mu     =  size(mu,1);
theta_a  =  zeros(n_mu, FOM.Qa);
theta_f  =  zeros(n_mu, FOM.Qf);
redMESH     =  FOM.HyRED.MESH;

for i = 1 : n_mu
    
    DATA_tmp            =  FOM.DATA;
    DATA_tmp.param      =  mu(i,:);
    DATA_tmp.Stiffening_power = mu(i,3);
    
    [~, A]              =  CSM_Assembler('internal_forces', redMESH, DATA_tmp, FOM.HyRED.FE_SPACE);
    [A_in, F_in]        =  CSM_ApplyBC(A, [], FOM.HyRED.FE_SPACE, redMESH, DATA_tmp);
    
    A_in  = A_in(:);
    
    % interpolation
    theta_a(i,1:Ma)     = ( FOM.HyRED.U_matrix(1:Ma, 1:Ma) \ A_in(FOM.HyRED.IDEIM_m(1:Ma))   ).';
    theta_f(i,1:Mf)     = ( FOM.HyRED.U_rhs(1:Mf, 1:Mf)    \ F_in(FOM.HyRED.IDEIM_rhs(1:Mf)) ).';
    
end

end

