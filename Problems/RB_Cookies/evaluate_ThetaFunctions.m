function [ theta_a, theta_f ] = evaluate_ThetaFunctions( mu, varargin )
  %   For reference, see Section 7.5 of
  %
  %   Quarteroni, Manzoni, Negri - REDUCED BASIS METHODS FOR PARTIAL
  %   DIFFERENTIAL EQUATIONS. AN INTRODUCTION. Springer, 2015

  %   This file is part of redbKIT.
  %   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
  %   Author: Federico Negri <federico.negri at epfl.ch>
  
Qa = 2;
Qf = 7;

n_mu     =  size(mu,1);
theta_a  =  zeros(n_mu, Qa);
theta_f  =  zeros(n_mu, Qf);

for i = 1 : n_mu

    theta_a(i,1) = 1;
    theta_a(i,2) = mu(i,1);

    theta_f(i,1) = mu(i,2);
    theta_f(i,2) = mu(i,3);
    theta_f(i,3) = mu(i,4);
    theta_f(i,4) = mu(i,5);
    theta_f(i,5) = mu(i,6);
    theta_f(i,6) = mu(i,7);
    theta_f(i,7) = mu(i,8);

end


end
