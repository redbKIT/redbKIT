function [ theta_a, theta_f ] = evaluate_ThetaFunctions( mu, FOM, Ma, Mf, nodes, ROM_SEMMT )
%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

if nargin < 3 || isempty(Ma)
    Ma = FOM.Qa;
end

if nargin < 4  || isempty(Mf)
    Mf = FOM.Qf;
end

n_mu     =  size(mu,1);
theta_a  =  zeros(n_mu, FOM.Qa);
theta_f  =  zeros(n_mu, FOM.Qf);

if nargin < 6  || isempty(ROM_SEMMT)
    loaded = load('../MeshMotion/ROM_SEMMT.mat');
    ROM_SEMMT = loaded.ROM_SEMMT; 
end

ROM_SEMMT.u_D{1} = @(x,mu)(deform_boundary(x(1,:),x(2,:),[],[1 mu]));
ROM_SEMMT.u_D{2} = @(x,mu)(deform_boundary(x(1,:),x(2,:),[],[2 mu]));

for i = 1 : n_mu
     
    MESH_tmp            =  FOM.HyRED.MESH;
    
    if nargin == 5
        MESH_tmp.nodes  = nodes;
    else
        MESH_tmp.nodes      =  MoveMesh(MESH_tmp.nodes, mu(i,FOM.DATA.shape_param), ROM_SEMMT);
    end
     
    MESH_tmp.vertices   =  MESH_tmp.nodes(:,1:MESH_tmp.numVertices);
    [MESH_tmp.jac, MESH_tmp.invjac, MESH_tmp.h] = geotrasf(MESH_tmp.dim, MESH_tmp.vertices, MESH_tmp.elements);   
    DATA_tmp            =  FOM.DATA;
    DATA_tmp.param      =  mu(i,:);
    
    [A, F]              =  ADR_Assembler(MESH_tmp, DATA_tmp, FOM.HyRED.FE_SPACE);
    [A_in, F_in]        =  ADR_ApplyBC(A, F, FOM.HyRED.FE_SPACE, MESH_tmp, DATA_tmp);
     
    A_in  = A_in(:);
     
    % interpolation
    theta_a(i,1:Ma)     = ( FOM.HyRED.U_matrix(1:Ma, 1:Ma) \ A_in(FOM.HyRED.IDEIM_m(1:Ma))   ).';
    theta_f(i,1:Mf)     = ( FOM.HyRED.U_rhs(1:Mf, 1:Mf)    \ F_in(FOM.HyRED.IDEIM_rhs(1:Mf)) ).';
    
end

end

