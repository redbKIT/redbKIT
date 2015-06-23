function [phi, dphix, dphiy] = fem_basis2D(fem, x, y, boundary)
%FEM_BASIS2D finite element basis functions.
%       F. Saleri 13-01-03,  F. Negri 22.11.2014


switch fem
    
    case 'P1' 
        
        phi(1,:) = 1-x-y;
        phi(2,:) = x;
        phi(3,:) = y;
        dphix(1,:) = -1+0.*x;
        dphix(2,:) =  1+0.*x;
        dphix(3,:) =  0+0.*x;
        dphiy(1,:) = -1+0.*x;
        dphiy(2,:) =  0+0.*x;
        dphiy(3,:) =  1+0.*x;  
        
    case 'P2' 
       
        phi(1,:) = (1-x-y).*(1-2*x-2*y);
        phi(2,:) = x.*(-1+2*x);
        phi(3,:) = y.*(-1+2*y);
        phi(4,:) = 4*x.*(1-x-y);
        phi(5,:) = 4*x.*y;
        phi(6,:) = 4*y.*(1-x-y);
        dphix(1,:) = -3+4*x+4*y;
        dphix(2,:) = -1+4*x;
        dphix(3,:) = 0 + 0.*x;
        dphix(4,:) = 4-8*x-4*y;
        dphix(5,:) = 4*y;
        dphix(6,:) = -4*y;
        dphiy(1,:) = -3+4*x+4*y;
        dphiy(2,:) = 0+0.*x;
        dphiy(3,:) = -1+4*y;
        dphiy(4,:) = -4*x;
        dphiy(5,:) = 4*x;
        dphiy(6,:) = 4-4*x-8*y;
    
    otherwise
        msg = [fem, ' finite element space not available'];
        error(msg)
end


if nargin == 4 && boundary == 1

    switch fem
        
        case 'P1'
            BDOF = [1 2];
        case 'P2'
            BDOF = [1 2 4];
        otherwise
            error(['Boundary integrals not supported for ', fem ,' FE space'])
    end
    phi   = phi(BDOF,:);
    dphix = dphix(BDOF,:);
    dphiy = dphiy(BDOF,:);
end

return
