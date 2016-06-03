function [ theta_a, theta_f ] = evaluate_ThetaFunctions( mu, FOM, M, varargin )

Qa = 3;

load EIM_quad EIM_quad;
load PHI_EIM PHI_EIM;

if nargin < 3
    M  = size(EIM_quad, 1); 
end

M        =  min(size(EIM_quad, 1), M);

Qf       =  size(EIM_quad, 1);

n_mu     =  size(mu,1);
theta_a  =  zeros(n_mu, Qa);
theta_f  =  zeros(n_mu, Qf);

for i = 1 : n_mu
    
    theta_a(i,1) = 1;
    theta_a(i,2) = cos(pi/180*(mu(i,3)));
    theta_a(i,3) = sin(pi/180*(mu(i,3)));  
    
    % evaluate source term in the quadrature nodes selected by EIM
    f_i_IEIM     = nonaffine_source(EIM_quad, mu);    
    
    % compute interpolation coefficient
    theta_f(i,1:M) = PHI_EIM(1:M,1:M) \ f_i_IEIM(1:M);
    
end


end

