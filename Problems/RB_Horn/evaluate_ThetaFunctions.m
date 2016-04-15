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

for i = 1 : n_mu
    
    MESH_tmp            =  FOM.HyRED.MESH;
    MESH_tmp.vertices   =  RBF_DeformGeometry(mu(i,FOM.DATA.shape_param),MESH_tmp.vertices);
    [MESH_tmp.jac, MESH_tmp.invjac, MESH_tmp.h] = geotrasf(MESH_tmp.dim, MESH_tmp.vertices, MESH_tmp.elements);   
    DATA_tmp            =  FOM.DATA;
    DATA_tmp.param      =  mu(i,:);
    
    [A, F]              =  Assembler_2D(MESH_tmp, DATA_tmp, FOM.HyRED.FE_SPACE);
    [A_in, F_in]        =  ADR_ApplyBC(A, F, FOM.HyRED.FE_SPACE, MESH_tmp, DATA_tmp);
    
    A_in  = A_in(:);
    
    
    % interpolation
    theta_a(i,1:Ma)     = ( FOM.HyRED.U_matrix(1:Ma, 1:Ma) \ A_in(FOM.HyRED.IDEIM_m(1:Ma))   ).';
    theta_f(i,1:Mf)     = ( FOM.HyRED.U_rhs(1:Mf, 1:Mf)    \ F_in(FOM.HyRED.IDEIM_rhs(1:Mf)) ).';
    
end

end

