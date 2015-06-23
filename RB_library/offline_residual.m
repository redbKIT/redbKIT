function [Cqq, dqq, Eqq] = offline_residual(FOM, V)
%OFFLINE_RESIDUAL offline computation of the mu-independent terms of the 
%dual norm of the residual 
%   
%   [CQQ, DQQ, EQQ] = OFFLINE_RESIDUAL(FOM, V) given a FOM struct and a
%   trial basis V, computes the mu-independent terms of the 
%   dual norm of the residual. CQQ is a ROM.Qf x ROM.Qf cell array of numbers,  
%   DQQ is a ROM.Qf x ROM.Qa cell array of ROM.N x 1 vectors,  while EQQ is 
%   a ROM.Qa x ROM.Qa cell array of ROM.N x ROM.N matrices.

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Fédérale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

for q1 = 1 : FOM.Qf
    
    t = FOM.Xnorm\FOM.Fq{q1};
    
    for q2 = 1: FOM.Qf
        
        Cqq{q1,q2} = t'*FOM.Fq{q2};
        
    end
    
end


for q1 = 1 : FOM.Qa
    
    Z = FOM.Xnorm\(FOM.Aq{q1}*V);
    
    for q2 = 1 : FOM.Qa
        Eqq{q1,q2} = Z'*(FOM.Aq{q2}*V);
    end
    
    for q2 = 1 : FOM.Qf
        dqq{q1,q2} = Z'*FOM.Fq{q2};
    end
end



end