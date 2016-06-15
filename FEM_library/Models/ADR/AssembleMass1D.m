function [M_restricted, M] = AssembleMass1D( fem, boundaries, vertices, nov, b_side, b_vertices)
%AssembleMass1D assemble 1D (as 2D boundary) mass

[~,nbn]       = select(fem, 2);

[csi,w]       = xwgl(6, 0, 1);
[phi]         = fem_basis(2, fem, [csi; 0*csi], 1);

one           = ones(nbn, 1);

M    = sparse(nov,nov);
MLOC = (phi.*w(one,:))*phi';

for side = b_side
    dof              = boundaries(1:nbn,side);
    x                = vertices(1,dof(1:2));
    y                = vertices(2,dof(1:2));
    side_length      = sqrt((x(2)-x(1))^2+(y(2)-y(1))^2);
    M(dof,dof)       = M(dof,dof) + side_length*MLOC; 
end

M_restricted = M(b_vertices, b_vertices);

return




