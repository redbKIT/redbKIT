function [ theta_a, theta_f ] = evaluate_ThetaFunctions( mu, varargin )

Qa = 4;
Qf = 1;

n_mu     =  size(mu,1);
theta_a  =  zeros(n_mu, Qa);
theta_f  =  zeros(n_mu, Qf);

for i = 1 : n_mu
    
    kappa = 2*pi*mu(i,1)/340;
    
    theta_a(i,1) = 1;
    theta_a(i,2) = -kappa^2;
    theta_a(i,3) = complex(0,1)*kappa;
    theta_a(i,4) = 1/2;
    
    theta_f(i,1) = 2*complex(0,1)*kappa;
    
end


end

