function [nln,nbn] = select2D(fem)
%SELECT2D number of local dof of a finite element 2D-space
%
% NLN is the number of local dofs on the reference triangle
% NBN is the number of local boundary dofs

switch fem
    
    case 'P1'
        nln = 3;
        nbn = 2;
    
    case 'P2'
        
        nln = 6;
        nbn = 3;
        
    otherwise
        
        msg = [fem, ' finite element space not available'];
        error(msg)
end

return