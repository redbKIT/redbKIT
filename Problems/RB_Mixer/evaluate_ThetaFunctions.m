function [ theta_a, theta_f ] = evaluate_ThetaFunctions( mu, varargin )
%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

Qa = 2;
Qf = 6;

n_mu     =  size(mu,1);
theta_a  =  zeros(n_mu, Qa);
theta_f  =  zeros(n_mu, Qf);

for i = 1 : n_mu
    
    theta_a(i,1) = 1/mu(i,4);
    theta_a(i,2) = 1;
    
    theta_f(i,1) = 1/mu(i,4)*mu(i,1);
    theta_f(i,2) = 1/mu(i,4)*mu(i,2);
    theta_f(i,3) = 1/mu(i,4)*mu(i,3);
    theta_f(i,4) = mu(i,1);
    theta_f(i,5) = mu(i,2);
    theta_f(i,6) = mu(i,3);
    
end


end