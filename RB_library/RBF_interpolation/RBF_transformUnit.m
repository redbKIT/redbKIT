function [mu_Unit]  =  RBF_transformUnit(mu,transform)
%RBF_TRANSFORMUNIT
%
%   [mu_Unit]  =  RBF_TRANSFORMUNIT(mu,transform)

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

P       = size(mu,2);
mu_Unit = 0*mu;

for i = 1 : P
      
      mu_Unit(:,i) = transform(i,1) + transform(i,2)*mu(:,i);
      
end



end