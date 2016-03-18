function [ deltaN ] = error_estimate(ROM, uN, mu, varargin)
%ERROR_ESTIMATE compute an estimate of the RB error for a given parameter
%value
%  
%   [ DELTAN ] = ERROR_ESTIMATE(ROM, UN, MU) given a ROM struct and the RB 
%   solution UN corresponding to the parameter value MU, computes an
%   estimate of the error between the reduced and full order solutions.

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

N = length(uN);

[ theta_a, theta_f ] = evaluate_ThetaFunctions( mu, ROM, varargin{:} );

INDX = 1:N;

%% evaluate a lower bound or an approximation to the stability factor
if isfield(ROM, 'stabFactor')
    beta   = RBF_OnlineInterpolation(ROM, mu);
else
    beta   = 1;
end

%% evaluate dual norm of the residual
res_aa = 0;
res_af = 0;
res_ff = 0;

for q1 = 1 : ROM.Qf
    for q2 = 1 : ROM.Qf
        res_ff = res_ff + theta_f(q1)' * theta_f(q2) * ROM.Cqq{q1,q2};
    end
end


for q1 = 1 : ROM.Qa
    for q2 = 1 : ROM.Qa
        res_aa = res_aa + theta_a(q1)' * theta_a(q2) * uN'*(ROM.Eqq{q1,q2}(INDX,INDX)*uN);
    end

    for q2 = 1 : ROM.Qf
        res_af = res_af + theta_a(q1)' * theta_f(q2)  * uN'*ROM.dqq{q1,q2}(INDX) ...
                        + theta_a(q1)  * theta_f(q2)' * ((ROM.dqq{q1,q2}(INDX))'*uN);
                        %+ theta_a(q1) * theta_f(q2)' * (ROM.dqq{2,q2,q1}(INDX)*uN);

    end
end


res = res_aa - (res_af) + res_ff;

deltaN = sqrt(abs(res)) / beta;

end