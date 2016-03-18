function [uh] = solve_HFsystem(FOM, mu)
%SOLVE_HFSYSTEM solves the high-fidelity (full order) problem
%
%   [UH] = SOLVE_HFSYSTEM(FOM, MU) given a parameter value MU, the function 
%   assembles the high-fidelity problem starting from the mu-independent 
%   arrays stored in FOM, and solves it.

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

[ theta_a, theta_f ] = evaluate_ThetaFunctions( mu, FOM );

Ah = theta_a(1) * FOM.Aq{1};
Fh = theta_f(1) * FOM.Fq{1};

for q = 2 : FOM.Qa
    Ah = Ah + theta_a(q) * FOM.Aq{q};
end

for q  = 2 : FOM.Qf
    Fh = Fh + theta_f(q) * FOM.Fq{q};
end

uh                            = zeros(FOM.MESH.numNodes,1);
uh(FOM.MESH.internal_dof)     = Ah \ Fh;
uh(FOM.MESH.Dirichlet_dof)    = FOM.u_D(FOM.MESH.nodes(:,FOM.MESH.Dirichlet_dof), mu);


end