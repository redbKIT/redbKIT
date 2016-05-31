function [rowdG, coldG, coefdG, rowG, coefG] = CSM_assembler_M_omp...
    (dim, Material_Model, material_param, Uh, elements, nln, w, invjac, ...
    detjac, phi, dphi_ref)%#codegen
%CSM_assembler_M_omp 2D/3D assembler for Hyper-Elasticity equations.

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri@epfl.ch>

noe          = size(elements,2);
nln2         = nln*nln;

n_blocks     = dim;
n_blocks2    = n_blocks^2;

rowdG         = zeros(noe*n_blocks2*nln2,1);
coldG         = rowdG;
coefdG        = rowdG;

rowG         = zeros(noe*n_blocks*nln,1);
coefG        = rowG;

NumQuadPoints          = length(w);

Id = eye(dim,dim);

NumNodes = length( Uh ) / dim;

global_index = 1;
global_indexG = 1;

switch Material_Model
    
    case 2 %St. Venant Kirchhoff
        
        Young = material_param(1);
        Poisson = material_param(2);
        mu = Young / (2 + 2 * Poisson);
        lambda =  Young * Poisson /( (1+Poisson) * (1-2*Poisson) );
        
        switch dim
            
            case 2
                for ie = 1 : noe
                    
                    e   = elements(1:nln,ie);
                    gradx = invjac(ie,1,1) * dphi_ref(:, :, 1) + invjac(ie,1,2) * dphi_ref(:, :, 2);
                    grady = invjac(ie,2,1) * dphi_ref(:, :, 1) + invjac(ie,2,2) * dphi_ref(:, :, 2);
                    
                    uhdx = Uh(e)' * gradx;
                    uhdy = Uh(e)' * grady;
                    
                    vhdx = Uh(e+NumNodes)' * gradx;
                    vhdy = Uh(e+NumNodes)' * grady;
                    
                    for i = 1 : n_blocks
                        
                        ei    = zeros(n_blocks,1);
                        ei(i) = 1;
                        eim   = repmat(ei,1,n_blocks);
                        
                        for j = 1 : n_blocks
                            
                            ej    = zeros(n_blocks,1);
                            ej(j) = 1;
                            ejm    = repmat(ej,1,n_blocks);
                            
                            for A = 1 : nln
                                
                                for B = 1 : nln
                                    
                                    intgr = 0;
                                    
                                    for q = 1 : NumQuadPoints
                                        
                                        GRAD_Uh = [uhdx(q) uhdy(q); vhdx(q) vhdy(q)];
                                        
                                        GRAD_V  = [gradx(A,q) grady(A,q);...
                                            gradx(A,q) grady(A,q)].*eim;
                                        
                                        GRAD_U  = [gradx(B,q) grady(B,q);...
                                            gradx(B,q) grady(B,q)].*ejm;
                                        
                                        F  = Id + GRAD_Uh;
                                        E  = 0.5 * (F'*F - Id);
                                        
                                        dF = GRAD_U;
                                        dE = 0.5 * (dF'  * F + F' * dF);
                                        
                                        dP =  dF * (2*mu*E  + lambda  * Trace(E)*Id) + ...
                                            F  * (2*mu*dE + lambda * Trace(dE)*Id);
                                        
                                        intgr  = intgr + ( Mdot( GRAD_V, dP ) ) * w(q);
                                        
                                    end
                                    
                                    coefdG(global_index) = intgr*detjac(ie);
                                    rowdG(global_index)  = elements(A,ie)+(i-1)*NumNodes;
                                    coldG(global_index)  = elements(B,ie)+(j-1)*NumNodes;
                                    global_index         = global_index + 1;
                                    
                                end
                            end
                        end
                    end
                    
                    %rhs
                    for i = 1 : n_blocks
                        
                        ei    = zeros(n_blocks,1);
                        ei(i) = 1;
                        eim   = repmat(ei,1,n_blocks);
                        
                        for A = 1 : nln
                            
                            intgr = 0;
                            
                            for q = 1 : NumQuadPoints
                                
                                GRAD_Uh = [uhdx(q) uhdy(q); vhdx(q) vhdy(q)];
                                
                                GRAD_V  = [gradx(A,q) grady(A,q);...
                                    gradx(A,q) grady(A,q)].*eim;
                                
                                F          =   Id + GRAD_Uh;
                                E          =   0.5 * (F'*F - Id);
                                P_Uh       =   F * ( 2*mu*E + lambda*Trace(E)*Id);
                                
                                intgr  = intgr + ( Mdot( GRAD_V, P_Uh ) ) * w(q);
                                
                            end
                            
                            coefG(global_indexG) = intgr*detjac(ie);
                            rowG(global_indexG)  = elements(A,ie)+(i-1)*NumNodes;
                            global_indexG        = global_indexG + 1;
                            
                        end
                        
                    end
                    
                end
                
            case 3
                
                for ie = 1 : noe
                    
                    e   = elements(1:nln,ie);
                    
                    gradx = invjac(ie,1,1) * dphi_ref(:, :, 1) + invjac(ie,1,2) * dphi_ref(:, :, 2) + invjac(ie,1,3) * dphi_ref(:, :, 3);
                    grady = invjac(ie,2,1) * dphi_ref(:, :, 1) + invjac(ie,2,2) * dphi_ref(:, :, 2) + invjac(ie,2,3) * dphi_ref(:, :, 3);
                    gradz = invjac(ie,3,1) * dphi_ref(:, :, 1) + invjac(ie,3,2) * dphi_ref(:, :, 2) + invjac(ie,3,3) * dphi_ref(:, :, 3);
                    
                    uhdx = Uh(e)' * gradx;
                    uhdy = Uh(e)' * grady;
                    uhdz = Uh(e)' * gradz;
                    
                    vhdx = Uh(e+NumNodes)' * gradx;
                    vhdy = Uh(e+NumNodes)' * grady;
                    vhdz = Uh(e+NumNodes)' * gradz;
                    
                    whdx = Uh(e+2*NumNodes)' * gradx;
                    whdy = Uh(e+2*NumNodes)' * grady;
                    whdz = Uh(e+2*NumNodes)' * gradz;
                    
                    
                    for i = 1 : n_blocks
                        
                        ei    = zeros(n_blocks,1);
                        ei(i) = 1;
                        eim   = repmat(ei,1,n_blocks);
                        
                        for j = 1 : n_blocks
                            
                            ej    = zeros(n_blocks,1);
                            ej(j) = 1;
                            ejm    = repmat(ej,1,n_blocks);
                            
                            for A = 1 : nln
                                
                                for B = 1 : nln
                                    
                                    intgr = 0;
                                    
                                    for q = 1 : NumQuadPoints
                                        
                                        GRAD_Uh = [uhdx(q) uhdy(q) uhdz(q); vhdx(q) vhdy(q) vhdz(q); whdx(q) whdy(q) whdz(q)];
                                        
                                        GRAD_V  = [gradx(A,q) grady(A,q) gradz(A,q);...
                                            gradx(A,q) grady(A,q) gradz(A,q);....
                                            gradx(A,q) grady(A,q) gradz(A,q)].*eim;
                                        
                                        GRAD_U  = [gradx(B,q) grady(B,q) gradz(B,q);...
                                            gradx(B,q) grady(B,q) gradz(B,q);....
                                            gradx(B,q) grady(B,q) gradz(B,q)].*ejm;
                                        
                                        F  = Id + GRAD_Uh;
                                        E  = 0.5 * (F'*F - Id);
                                        
                                        dF = GRAD_U;
                                        dE = 0.5 * (dF'  * F + F' * dF);
                                        
                                        dP =  dF * (2*mu*E  + lambda  * Trace(E)*Id) + ...
                                            F  * (2*mu*dE + lambda * Trace(dE)*Id);
                                        
                                        intgr  = intgr + ( Mdot( GRAD_V, dP ) ) * w(q);
                                        
                                    end
                                    
                                    coefdG(global_index) = intgr*detjac(ie);
                                    rowdG(global_index)  = elements(A,ie)+(i-1)*NumNodes;
                                    coldG(global_index)  = elements(B,ie)+(j-1)*NumNodes;
                                    global_index         = global_index + 1;
                                    
                                end
                            end
                        end
                    end
                    
                    %rhs
                    for i = 1 : n_blocks
                        
                        ei    = zeros(n_blocks,1);
                        ei(i) = 1;
                        eim   = repmat(ei,1,n_blocks);
                        
                        for A = 1 : nln
                            
                            intgr = 0;
                            
                            for q = 1 : NumQuadPoints
                                
                                GRAD_Uh = [uhdx(q) uhdy(q) uhdz(q); vhdx(q) vhdy(q) vhdz(q); whdx(q) whdy(q) whdz(q)];
                                
                                GRAD_V  = [gradx(A,q) grady(A,q) gradz(A,q);...
                                    gradx(A,q) grady(A,q) gradz(A,q);....
                                    gradx(A,q) grady(A,q) gradz(A,q)].*eim;
                                
                                F          =   Id + GRAD_Uh;
                                E          =   0.5 * (F'*F - Id);
                                P_Uh       =   F * ( 2*mu*E + lambda*Trace(E)*Id);
                                
                                intgr  = intgr + ( Mdot( GRAD_V, P_Uh ) ) * w(q);
                                
                            end
                            
                            coefG(global_indexG) = intgr*detjac(ie);
                            rowG(global_indexG)  = elements(A,ie)+(i-1)*NumNodes;
                            global_indexG        = global_indexG + 1;
                            
                        end
                        
                    end
                    
                end
        end
        
        
    case 3 %Neo-Hookean
        
        Young = material_param(1);
        Poisson = material_param(2);
        mu = Young / (2 + 2 * Poisson);
        lambda =  Young * Poisson /( (1+Poisson) * (1-2*Poisson) );
        
        switch dim
            
            case 2
                for ie = 1 : noe
                    
                    e   = elements(1:nln,ie);
                    gradx = invjac(ie,1,1) * dphi_ref(:, :, 1) + invjac(ie,1,2) * dphi_ref(:, :, 2);
                    grady = invjac(ie,2,1) * dphi_ref(:, :, 1) + invjac(ie,2,2) * dphi_ref(:, :, 2);
                    
                    uhdx = Uh(e)' * gradx;
                    uhdy = Uh(e)' * grady;
                    
                    vhdx = Uh(e+NumNodes)' * gradx;
                    vhdy = Uh(e+NumNodes)' * grady;
                    
                    for i = 1 : n_blocks
                        
                        ei    = zeros(n_blocks,1);
                        ei(i) = 1;
                        eim   = repmat(ei,1,n_blocks);
                        
                        for j = 1 : n_blocks
                            
                            ej    = zeros(n_blocks,1);
                            ej(j) = 1;
                            ejm    = repmat(ej,1,n_blocks);
                            
                            for A = 1 : nln
                                
                                for B = 1 : nln
                                    
                                    intgr = 0;
                                    
                                    for q = 1 : NumQuadPoints
                                        
                                        GRAD_Uh = [uhdx(q) uhdy(q); vhdx(q) vhdy(q)];
                                        
                                        GRAD_V  = [gradx(A,q) grady(A,q);...
                                            gradx(A,q) grady(A,q)].*eim;
                                        
                                        GRAD_U  = [gradx(B,q) grady(B,q);...
                                            gradx(B,q) grady(B,q)].*ejm;
                                        
                                        F  = Id + GRAD_Uh;
                                        invFT      =   inv(F)';
                                        J          =   abs(det(F));
                                        
                                        dF         = GRAD_U;
                                        
                                        dP = mu * dF + (mu - lambda * log(J)) * invFT * dF.' * invFT ...
                                            + lambda * Trace(invFT.' * dF)*invFT;
                                        
                                        intgr  = intgr + ( Mdot( GRAD_V, dP ) ) * w(q);
                                        
                                    end
                                    
                                    coefdG(global_index) = intgr*detjac(ie);
                                    rowdG(global_index)  = elements(A,ie)+(i-1)*NumNodes;
                                    coldG(global_index)  = elements(B,ie)+(j-1)*NumNodes;
                                    global_index         = global_index + 1;
                                    
                                end
                            end
                        end
                    end
                    
                    %rhs
                    for i = 1 : n_blocks
                        
                        ei    = zeros(n_blocks,1);
                        ei(i) = 1;
                        eim   = repmat(ei,1,n_blocks);
                        
                        for A = 1 : nln
                            
                            intgr = 0;
                            
                            for q = 1 : NumQuadPoints
                                
                                GRAD_Uh = [uhdx(q) uhdy(q); vhdx(q) vhdy(q)];
                                
                                GRAD_V  = [gradx(A,q) grady(A,q);...
                                    gradx(A,q) grady(A,q)].*eim;
                                
                                F          =   Id + GRAD_Uh;
                                invFT      =   inv(F)';
                                J          =   abs(det(F));
                                
                                P_Uh       =   mu * (F - invFT) + lambda * log(J) * invFT;
                                
                                intgr  = intgr + ( Mdot( GRAD_V, P_Uh ) ) * w(q);
                                
                            end
                            
                            coefG(global_indexG) = intgr*detjac(ie);
                            rowG(global_indexG)  = elements(A,ie)+(i-1)*NumNodes;
                            global_indexG        = global_indexG + 1;
                            
                        end
                        
                    end
                    
                end
                
            case 3
                
                for ie = 1 : noe
                    
                    e   = elements(1:nln,ie);
                    
                    gradx = invjac(ie,1,1) * dphi_ref(:, :, 1) + invjac(ie,1,2) * dphi_ref(:, :, 2) + invjac(ie,1,3) * dphi_ref(:, :, 3);
                    grady = invjac(ie,2,1) * dphi_ref(:, :, 1) + invjac(ie,2,2) * dphi_ref(:, :, 2) + invjac(ie,2,3) * dphi_ref(:, :, 3);
                    gradz = invjac(ie,3,1) * dphi_ref(:, :, 1) + invjac(ie,3,2) * dphi_ref(:, :, 2) + invjac(ie,3,3) * dphi_ref(:, :, 3);
                    
                    uhdx = Uh(e)' * gradx;
                    uhdy = Uh(e)' * grady;
                    uhdz = Uh(e)' * gradz;
                    
                    vhdx = Uh(e+NumNodes)' * gradx;
                    vhdy = Uh(e+NumNodes)' * grady;
                    vhdz = Uh(e+NumNodes)' * gradz;
                    
                    whdx = Uh(e+2*NumNodes)' * gradx;
                    whdy = Uh(e+2*NumNodes)' * grady;
                    whdz = Uh(e+2*NumNodes)' * gradz;
                    
                    
                    for i = 1 : n_blocks
                        
                        ei    = zeros(n_blocks,1);
                        ei(i) = 1;
                        eim   = repmat(ei,1,n_blocks);
                        
                        for j = 1 : n_blocks
                            
                            ej    = zeros(n_blocks,1);
                            ej(j) = 1;
                            ejm    = repmat(ej,1,n_blocks);
                            
                            for A = 1 : nln
                                
                                for B = 1 : nln
                                    
                                    intgr = 0;
                                    
                                    for q = 1 : NumQuadPoints
                                        
                                        GRAD_Uh = [uhdx(q) uhdy(q) uhdz(q); vhdx(q) vhdy(q) vhdz(q); whdx(q) whdy(q) whdz(q)];
                                        
                                        GRAD_V  = [gradx(A,q) grady(A,q) gradz(A,q);...
                                            gradx(A,q) grady(A,q) gradz(A,q);....
                                            gradx(A,q) grady(A,q) gradz(A,q)].*eim;
                                        
                                        GRAD_U  = [gradx(B,q) grady(B,q) gradz(B,q);...
                                            gradx(B,q) grady(B,q) gradz(B,q);....
                                            gradx(B,q) grady(B,q) gradz(B,q)].*ejm;
                                        
                                        F  = Id + GRAD_Uh;
                                        invFT      =   inv(F)';
                                        J          =   abs(det(F));
                                        
                                        dF         = GRAD_U;
                                        
                                        dP = mu * dF + (mu - lambda * log(J)) * invFT * dF.' * invFT ...
                                            + lambda * Trace(invFT.' * dF)*invFT;
                                        
                                        intgr  = intgr + ( Mdot( GRAD_V, dP ) ) * w(q);
                                        
                                    end
                                    
                                    coefdG(global_index) = intgr*detjac(ie);
                                    rowdG(global_index)  = elements(A,ie)+(i-1)*NumNodes;
                                    coldG(global_index)  = elements(B,ie)+(j-1)*NumNodes;
                                    global_index         = global_index + 1;
                                    
                                end
                            end
                        end
                    end
                    
                    %rhs
                    for i = 1 : n_blocks
                        
                        ei    = zeros(n_blocks,1);
                        ei(i) = 1;
                        eim   = repmat(ei,1,n_blocks);
                        
                        for A = 1 : nln
                            
                            intgr = 0;
                            
                            for q = 1 : NumQuadPoints
                                
                                GRAD_Uh = [uhdx(q) uhdy(q) uhdz(q); vhdx(q) vhdy(q) vhdz(q); whdx(q) whdy(q) whdz(q)];
                                
                                GRAD_V  = [gradx(A,q) grady(A,q) gradz(A,q);...
                                    gradx(A,q) grady(A,q) gradz(A,q);....
                                    gradx(A,q) grady(A,q) gradz(A,q)].*eim;
                                
                                F          =   Id + GRAD_Uh;
                                invFT      =   inv(F)';
                                J          =   abs(det(F));
                                
                                P_Uh       =   mu * (F - invFT) + lambda * log(J) * invFT;
                                
                                intgr  = intgr + ( Mdot( GRAD_V, P_Uh ) ) * w(q);
                                
                            end
                            
                            coefG(global_indexG) = intgr*detjac(ie);
                            rowG(global_indexG)  = elements(A,ie)+(i-1)*NumNodes;
                            global_indexG        = global_indexG + 1;
                            
                        end
                        
                    end
                    
                end
        end
end
end

%==========================================================================
function P = Mdot(X,Y)
P = 0;
for i = 1 : size(X, 1)
    for j = 1 : size(X, 1)
        P = P + X(i,j) * Y(i,j);
    end
end

end
%==========================================================================
function T = Trace(X)

T = 0;
for i = 1 : size(X,2)
    T = T + X(i,i);
end

end
%==========================================================================
function [gradx, grady] = Gradient_BasisFunctions2D(invjac, GradRef, ie)

gradx = invjac(ie,1,1) * GradRef(:, :, 1) + invjac(ie,1,2) * GradRef(:, :, 2);
grady = invjac(ie,2,1) * GradRef(:, :, 1) + invjac(ie,2,2) * GradRef(:, :, 2);

end
%==========================================================================
function [gradx, grady, gradz] = Gradient_BasisFunctions3D(invjac, GradRef, ie)

gradx = invjac(ie,1,1) * GradRef(:, :, 1) + invjac(ie,1,2) * GradRef(:, :, 2) + invjac(ie,1,3) * GradRef(:, :, 3);
grady = invjac(ie,2,1) * GradRef(:, :, 1) + invjac(ie,2,2) * GradRef(:, :, 2) + invjac(ie,2,3) * GradRef(:, :, 3);
gradz = invjac(ie,3,1) * GradRef(:, :, 1) + invjac(ie,3,2) * GradRef(:, :, 2) + invjac(ie,3,3) * GradRef(:, :, 3);

end
%==========================================================================
function [dP] = SVK_SigmaTensor(GRAD_Uh, GRAD_U, Id, mu, lambda)

F  = Id + GRAD_Uh;
E  = 0.5 * (F'*F - Id);

dF = GRAD_U;
dE = 0.5 * (dF'  * F + F' * dF);

dP =  dF * (2*mu*E  + lambda  * Trace(E)*Id) + F  * (2*mu*dE + lambda * Trace(dE)*Id);

end
%==========================================================================


%% Old STUFF

%
%         localcoefdG = zeros(1,n_blocks2*nln2);
%         localrowdG = zeros(1,n_blocks2*nln2);
%         localcoldG = zeros(1,n_blocks2*nln2);
%         localcoefG = zeros(1,n_blocks*nln);
%         localrowG  = zeros(1,n_blocks*nln);
%
%         gradphi = zeros(dim, nln, NumQuadPoints);
%         F       = zeros(dim, dim, NumQuadPoints);
%         E       = zeros(dim, dim, NumQuadPoints);
%         traceE  = zeros(1, NumQuadPoints);
%         GradUh  = zeros(dim, dim, NumQuadPoints);
%
%         for ie = 1 : noe
%
%             for q = 1 : NumQuadPoints
%                 for k = 1 : nln
%                     for d1 = 1 : dim
%                         gradphi(d1,k,q) = 0;
%                         for d2 = 1 : dim
%                             gradphi(d1,k,q) = gradphi(d1,k,q) + invjac(ie,d1,d2)*dphi_ref(k, q, d2);% check
%                         end
%                     end
%                 end
%
%                 for d1 = 1 : dim
%                     for d2 = 1 : dim
%                         for k = 1 : nln
%                             e_k = elements( k, ie ) + (d1-1)*NumNodes;
%                             GradUh(d1, d2, q) = GradUh(d1, d2, q) + Uh(e_k) * gradphi(d2, k, q);
%                         end
%                     end
%                 end
%
%                 F(:, :, q)  = Id + GradUh(:, :, q);
%                 E(:, :, q)  = 0.5 * (F(:,:,q)'*F(:,:,q) - Id);
%                 traceE(q)   = Trace( E(:,:,q) );
%             end
%
%
%             iii = 1;
%             ii  = 1;
%
%             for A = 1 : nln
%
%                 for i = 1 : n_blocks
%
%                     GRAD_V = zeros(dim, dim);
%
%                     for B = 1 : nln
%
%                         for j = 1 : n_blocks
%
%                             GRAD_U = zeros(dim, dim);
%
%                             intgr = 0;
%
%                             for q = 1 : NumQuadPoints
%
%                                 for d2 = 1 : dim
%                                     GRAD_V(i, d2) = gradphi(d2, A, q);
%                                     GRAD_U(j, d2) = gradphi(d2, B, q);
%                                 end
%
%                                 dF = GRAD_U;
%                                 dE = 0.5 * (dF'  * F(:,:,q) + F(:,:,q)' * dF);
%
%                                 dP =  dF * (2*mu*E(:,:,q)  + lambda  * traceE(q)*Id) + ...
%                                     F(:,:,q)  * (2*mu*dE + lambda * Trace(dE)*Id);
%
%                                 intgr  = intgr + ( Mdot( GRAD_V, dP ) ) * w(q);
%
%                             end
%
%                             localcoefdG(iii) = intgr*detjac(ie);
%                             localrowdG(iii)  = elements(A,ie)+(i-1)*NumNodes;
%                             localcoldG(iii)  = elements(B,ie)+(j-1)*NumNodes;
%                             iii         = iii + 1;
%                         end
%                     end
%
%                     intgr = 0;
%                     for q = 1 : NumQuadPoints
%
%                         for d2 = 1 : dim
%                              GRAD_V(i, d2) = gradphi(d2, A, q);
%                         end
%                         P_Uh    =   F(:,:,q) * ( 2*mu*E(:,:,q) + lambda*traceE(q)*Id);
%
%                         intgr  = intgr + ( Mdot( GRAD_V, P_Uh ) ) * w(q);
%                     end
%
%                     localcoefG(ii) = intgr*detjac(ie);
%                     localrowG(ii)  = elements(A,ie)+(i-1)*NumNodes;
%                     ii        = ii + 1;
%
%                 end
%             end
%
%             coefdG(ie, :)   = localcoefdG;
%             rowdG(ie, :)    = localrowdG;
%             coldG(ie, :)    = localcoldG;
%
%             coefG(ie, :)   = localcoefG;
%             rowG(ie, :)    = localrowG;
%         end
% end
%
% coefdG = coefdG(:);
% rowdG  = rowdG(:);
% coldG  = coldG(:);
%
% coefG  = coefG(:);
% rowG   = rowG(:);