function [def_nodes] = MoveMesh(nodes, mu, ROM_SEMMT, filename)

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Fédérale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

redMESH  =  ROM_SEMMT.HyRED.MESH;

DATA_tmp            =  ROM_SEMMT.DATA;
DATA_tmp.param      =  mu(1:3);
DATA_tmp.Stiffening_power = mu(3);

% [~, A]              =  CSM_Assembler_function('internal_forces', redMESH, DATA_tmp, ROM_SEMMT.HyRED.FE_SPACE);
% [A_in, F_in]        =  CSM_ApplyBC(A, [], ROM_SEMMT.HyRED.FE_SPACE, redMESH, DATA_tmp);

SolidModel          =  CSM_Assembler( redMESH, DATA_tmp, ROM_SEMMT.HyRED.FE_SPACE );
A                   =  SolidModel.compute_jacobian();
[A_in, F_in]        =  CSM_ApplyBC(A, [], ROM_SEMMT.HyRED.FE_SPACE, redMESH, DATA_tmp);

A_in  = A_in(:);

% interpolation
theta_a     = ( ROM_SEMMT.HyRED.U_matrix  \ A_in(ROM_SEMMT.HyRED.IDEIM_m)   ).';
theta_f     = ( ROM_SEMMT.HyRED.U_rhs     \ F_in(ROM_SEMMT.HyRED.IDEIM_rhs) ).';


AN = 0 * ROM_SEMMT.ANq{1};
FN = 0 * ROM_SEMMT.FNq{1};

for q = 1 : ROM_SEMMT.Qa
    AN = AN + theta_a(q) * ROM_SEMMT.ANq{q};
end

for q  = 1 : ROM_SEMMT.Qf
    FN = FN + theta_f(q) * ROM_SEMMT.FNq{q};
end

dN   = AN \ FN;

dNh                              = zeros(ROM_SEMMT.MESH.numNodes,1);
dNh(ROM_SEMMT.MESH.internal_dof) = ROM_SEMMT.V*dN;

for k = 1 : ROM_SEMMT.MESH.dim
    dNh(ROM_SEMMT.MESH.Dirichlet_dof_c{k}+(k-1)*ROM_SEMMT.MESH.numNodes)    = ROM_SEMMT.u_D{k}(ROM_SEMMT.MESH.nodes(:,ROM_SEMMT.MESH.Dirichlet_dof_c{k}), mu);
end
  
def_nodes = nodes+[dNh(1:ROM_SEMMT.MESH.numNodes)';dNh(1+ROM_SEMMT.MESH.numNodes:end)'];
   
if nargin == 4
    CSM_export_solution(2, dNh, ROM_SEMMT.MESH.vertices,...
        ROM_SEMMT.MESH.elements, ROM_SEMMT.MESH.numNodes, filename{1}, filename{2});
end

end