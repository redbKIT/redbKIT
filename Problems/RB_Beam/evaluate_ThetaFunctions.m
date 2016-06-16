function [ theta_a, theta_f ] = evaluate_ThetaFunctions( mu, varargin )
%   For reference, see Section 8.6 and 9.2 of
%
%   Quarteroni, Manzoni, Negri - REDUCED BASIS METHODS FOR PARTIAL
%   DIFFERENTIAL EQUATIONS. AN INTRODUCTION. Springer, 2015

Qa = 2;
Qf = 3;

n_mu     =  size(mu,1);
theta_a  =  zeros(n_mu, Qa);
theta_f  =  zeros(n_mu, Qf);


for i = 1 : n_mu
    
    theta_a(i,1) = mu(i,1)/(2+2*mu(i,2));
    theta_a(i,2) = mu(i,1)*mu(i,2)/((1+mu(i,2))*(1-2*mu(i,2)));
    
    theta_f(i,1) = 1e-4*mu(i,3);
    theta_f(i,2) = 1e-4*mu(i,4);
    theta_f(i,3) = 1e-4*mu(i,5);
    
end