function [d_F] = FSI_harmonicExtension(MESH, Displacement, HE_matrix)


d_F = zeros(MESH.Fluid.numNodes * MESH.dim, 1);

internal_dofs = HE_matrix.internal_dofs;

interface_dofs  = [];
d_S             = [];
Interface_SFmap = [];

for k = 1 : MESH.dim
    interface_dofs    = [interface_dofs; MESH.Fluid.numNodes*(k-1)+MESH.Fluid.dof_interface{k}];
    tmp               = Displacement(MESH.Solid.numNodes*(k-1)+MESH.Solid.dof_interface{k});
    d_S               = [d_S; tmp(MESH.Interface_SFmap{k})];
end


F   = - HE_matrix.HE(internal_dofs,interface_dofs) * d_S;

x   = HE_matrix.U \ (HE_matrix.L \ F(HE_matrix.perm));

d_F(internal_dofs) = x(HE_matrix.invp);

for k = 1 : MESH.dim
    
    tmp     = Displacement(MESH.Solid.numNodes*(k-1)+MESH.Solid.dof_interface{k});
    d_F(MESH.Fluid.numNodes*(k-1)+MESH.Fluid.dof_interface{k})  = tmp(MESH.Interface_SFmap{k});
    
end

d_F = reshape(d_F, MESH.Fluid.numNodes,  MESH.dim)';

end