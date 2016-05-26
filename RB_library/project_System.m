function [ANq, FNq] = project_System(FOM, V, method)
%PROJECT_SYSTEM project affine full order model over reduced bases to
%generate the reduced order model
%   
%   [AN, FN] = PROJECT_SYSTEM(FOM, V, METHOD) given a FOM struct containing
%   fields Aq, Fq, Qa and Qf, and a trial basis V, returns the reduced 
%   matrices ANq and FNq by either Galerkin (method = 'Galerkin') or Least-
%   Squares (method = 'LeastSquares') projection.

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

if nargin < 3 || isempty(method)
    method = 'Galerkin';
end


switch method
    case 'Galerkin'
        
        parfor q = 1 : FOM.Qa
            ANq{q} = V'*(FOM.Aq{q}*V);
        end
        
        parfor q = 1 : FOM.Qf
            FNq{q} = V'*FOM.Fq{q};
        end
        
    case 'LeastSquares'
        
        for q1 = 1 : FOM.Qa
            Z = FOM.Xnorm\(FOM.Aq{q1}*V);
            for q2 = 1 : FOM.Qa
                ANq{q1,q2} = Z'*(FOM.Aq{q2}*V);
            end
            
            for q2 = 1 : FOM.Qf
                FNq{q1,q2} = Z'*FOM.Fq{q2};
            end
        end
        
end