function [] = reset_randomSeed()
%RESET_RANDOMSEED wrapper to RANDSTREAM and RNG

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 


if verLessThan('matlab', '7.12')
    
    %type      = RandStream.DefaultStartupType;
    %s         = RandStream(type,'Seed',0);
    %RandStream.setGlobalStream(s);
    
    reset(RandStream.getDefaultStream);
    
else
    
    rng('default');
    
end


return