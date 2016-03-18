function [phi, dphi] = fem_basis(dim, fem, X, boundary)
%FEM_BASIS2D finite element basis functions.
%       F. Saleri 13-01-03,  F. Negri 22.11.2014

phi   = [];
dphix = [];
dphiy = [];
dphiz = [];

switch dim
    
    case 2
        
        x = X(1,:);
        y = X(2,:);
        
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
        
        dphi(:,:,1) = dphix;
        dphi(:,:,2) = dphiy;
        
    case 3
        
        x = X(1,:);
        y = X(2,:);
        z = X(3,:);
        
        switch fem
            case 'P1'
                phi(1,:) = 1-x-y-z;
                phi(2,:) = x;
                phi(3,:) = y;
                phi(4,:) = z;
                dphix(1,:) = -1+0.*x;
                dphix(2,:) =  1+0.*x;
                dphix(3,:) =  0+0.*x;
                dphix(4,:) =  0+0.*x;
                dphiy(1,:) = -1+0.*x;
                dphiy(2,:) =  0+0.*x;
                dphiy(3,:) =  1+0.*x;
                dphiy(4,:) =  0+0.*x;
                dphiz(1,:) = -1+0.*x;
                dphiz(2,:) =  0+0.*x;
                dphiz(3,:) =  0+0.*x;
                dphiz(4,:) =  1+0.*x;
            case 'P2'
                phi(1,:) = (1-x-y-z).*(1-2*x-2*y-2*z);
                phi(2,:) = x.*(2*x-1);
                phi(3,:) = y.*(2*y-1);
                phi(4,:) = z.*(2*z-1);
                phi(5,:) = 4*x.*(1-x-y-z);
                phi(6,:) = 4*x.*y;
                phi(7,:) = 4*y.*(1-x-y-z);
                phi(8,:) = 4*z.*(1-x-y-z);
                phi(9,:) = 4*x.*z;
                phi(10,:) = 4*y.*z;
                dphix(1,:) = -3+4*x+4*y+4*z;
                dphix(2,:) =  4*x-1;
                dphix(3,:) =  0+0.*x;
                dphix(4,:) =  0+0.*x;
                dphix(5,:) = 4-8*x-4*y-4*z;
                dphix(6,:) =  4*y;
                dphix(7,:) =  -4*y;
                dphix(8,:) =  -4*z;
                dphix(9,:) =  4*z;
                dphix(10,:) =  0+0.*x;
                dphiy(1,:) = -3+4*x+4*y+4*z;
                dphiy(2,:) =  0+0.*x;
                dphiy(3,:) =  -1+4*y;
                dphiy(4,:) =  0+0.*x;
                dphiy(5,:) = -4*x;
                dphiy(6,:) =  4*x;
                dphiy(7,:) =  4-4*x-8*y-4*z;
                dphiy(8,:) =  -4.*z;
                dphiy(9,:) =  0+0.*x;
                dphiy(10,:) =  4*z;
                dphiz(1,:) = -3+4*x+4*y+4*z;
                dphiz(2,:) =  0+0.*x;
                dphiz(3,:) =  0+0.*x;
                dphiz(4,:) =  -1+4*z;
                dphiz(5,:) = -4*x;
                dphiz(6,:) =  0+0.*x;
                dphiz(7,:) =  -4*y;
                dphiz(8,:) =  4-4*x-4*y-8*z;
                dphiz(9,:) =  4*x;
                dphiz(10,:) =  4*y;
                
            case 'B1' % FN
                phi(1,:) = (1-x-y-z).*(1-64*x.*y.*z);
                phi(2,:) = x.*(1-64.*y.*z.*(1-x-y-z));
                phi(3,:) = y.*(1-64.*x.*z.*(1-x-y-z));
                phi(4,:) = z.*(1-64.*x.*y.*(1-x-y-z));
                phi(5,:) = 256*x.*y.*z.*(1-x-y-z); % (d+1^(d+1) lambda_1 * lambda_2 * lambda_3)
                
                dphix(1,:) =  64.*y.*z.*(x + y + z - 1) + 64.*x.*y.*z - 1;
                dphix(2,:) =  64.*y.*z.*(x + y + z - 1) + 64.*x.*y.*z + 1;
                dphix(3,:) =  y.*(64.*x.*z + 64.*z.*(x + y + z - 1));
                dphix(4,:) =  z.*(64.*x.*y + 64.*y.*(x + y + z - 1));
                dphix(5,:) =  - 256.*y.*z.*(x + y + z - 1) - 256.*x.*y.*z;
                
                dphiy(1,:) =  64.*x.*z.*(x + y + z - 1) + 64.*x.*y.*z - 1;
                dphiy(2,:) =  x.*(64.*y.*z + 64.*z.*(x + y + z - 1));
                dphiy(3,:) =  64.*x.*z.*(x + y + z - 1) + 64.*x.*y.*z + 1;
                dphiy(4,:) =  z.*(64.*x.*y + 64.*x.*(x + y + z - 1));
                dphiy(5,:) =  - 256.*x.*z.*(x + y + z - 1) - 256.*x.*y.*z;
                
                dphiz(1,:) =  64.*x.*y.*(x + y + z - 1) + 64.*x.*y.*z - 1;
                dphiz(2,:) =  x.*(64.*y.*z + 64.*y.*(x + y + z - 1));
                dphiz(3,:) =  y.*(64.*x.*z + 64.*x.*(x + y + z - 1));
                dphiz(4,:) =  64.*x.*y.*(x + y + z - 1) + 64.*x.*y.*z + 1;
                dphiz(5,:) =  - 256.*x.*y.*(x + y + z - 1) - 256.*x.*y.*z;
                
            otherwise
                error([fem,' FE space is not available in 3D'])
                phi = []; dphix=[]; dphiy=[];
        end
        
        if nargin == 4 && boundary == 1
            
            switch fem
                case 'P1'
                    BDOF = [1 2 3];
                case 'P2'
                    BDOF = [1 2 3 5 6 7];
                case 'B1'
                    BDOF = [1 2 3];
                otherwise
                    error(['Boundary integrals not supported for ', fem ,' FE space'])
            end
            phi   = phi(BDOF,:);
            dphix = dphix(BDOF,:);
            dphiy = dphiy(BDOF,:);
            dphiz = dphiz(BDOF,:);
        end
        
        dphi(:,:,1) = dphix;
        dphi(:,:,2) = dphiy;
        dphi(:,:,3) = dphiz;
end

return
