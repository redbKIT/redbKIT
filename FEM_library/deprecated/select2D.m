function [nln,nbn] = select2D(fem)
%SELECT2D number of local dof of a finite element 2D-space
%
% NLN is the number of local dofs on the reference triangle
% NBN is the number of local boundary dofs

warning('select2D is deprecated. Use select instead.')

[nln,nbn] = select(fem, 2);

return
