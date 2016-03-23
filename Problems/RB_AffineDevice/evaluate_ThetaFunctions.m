function [ theta_a, theta_f ] = evaluate_ThetaFunctions( mu, varargin )

Qa = 6;
Qf = 1;

n_mu     =  size(mu,1);
theta_a  =  zeros(n_mu, Qa);
theta_f  =  zeros(n_mu, Qf);

for i = 1 : n_mu
    
    theta_a(i,1) = 1/(1+3*mu(i,1));
    theta_a(i,2) = 1 + 3*mu(i,1);
    theta_a(i,3) = mu(i,2)*1/((3*mu(i,1) + 1)^3/162) * (1 + 3*mu(i,1)).^3;
    theta_a(i,4) = 100;
    theta_a(i,5) = mu(i,3);
    theta_a(i,6) = 1;
    
    theta_f(i,1) = 10;
    
end


end

