function [errorL2,errorH1] = FEM_error(uh, MESH, DATA, FE_SPACE, t)
%FEM_ERROR compute L2, H1 and RHO error for 2D/3D fem approximation
%
%   [ERRORL2] = FEM_ERROR(UH, MESH, DATA, FE_SPACE) computes
%   the L^2 error for FEM finite elements solution UH. The exact solution
%   must be precised in the DATA.uexact field (see datafile template).
%
%   [ERRORL2,ERRORH1] = FEM_ERROR(UH, MESH, DATA, FE_SPACE) computes
%   the L^2 and the H^1 error. In this case the struct DATA must provide
%   also the field uxexact and uyexact for the first derivatives
%   of the exact solution.

%   Auhtors: F. Saleri 09-01-05. F. Negri 18.11.2014
%   3D by F. Negri

if nargin == 4
    t = [];
end

errorL2  = 0;
errorH1  = 0;

switch MESH.dim
    
    case 2
        
        x         = zeros(MESH.numElem,FE_SPACE.numQuadNodes);
        y         = x;
        coord_ref = MESH.chi;
        
        for j = 1:3
            i = MESH.elements(j,:);
            vtemp = MESH.vertices(1,i);
            x = x + vtemp'*coord_ref(j,:);
            vtemp = MESH.vertices(2,i);
            y = y + vtemp'*coord_ref(j,:);
        end
        
        
        uq      = DATA.uexact(x,y,t,DATA.param);
        
        
        uhe = uh(MESH.elements(1:FE_SPACE.numElemDof,:));
        uhq = uhe'*FE_SPACE.phi;
        
        if nargout == 1  % only L2 error is computed
            
            
            errorL2  = sum((((uhq-uq).^2)*FE_SPACE.quad_weights').*MESH.jac');
            
        elseif nargout == 2 % compute L2 and H1 errors
            
            uxq      = DATA.uxexact(x,y,t,DATA.param);
            uyq      = DATA.uyexact(x,y,t,DATA.param);
            
            dcdx     = MESH.invjac(:,1,1)';
            dedx     = MESH.invjac(:,1,2)';
            dcdy     = MESH.invjac(:,2,1)';
            dedy     = MESH.invjac(:,2,2)';
            
            
            dcdx     = dcdx(ones(FE_SPACE.numElemDof, 1), :);
            dedx     = dedx(ones(FE_SPACE.numElemDof, 1), :);
            uhdxq    = (uhe.*dcdx)'*FE_SPACE.dphi_ref(:,:,1)+(uhe.*dedx)'*FE_SPACE.dphi_ref(:,:,2);
            dcdy     = dcdy(ones(FE_SPACE.numElemDof, 1), :);
            dedy     = dedy(ones(FE_SPACE.numElemDof, 1), :);
            uhdyq    = (uhe.*dcdy)'*FE_SPACE.dphi_ref(:,:,1)+(uhe.*dedy)'*FE_SPACE.dphi_ref(:,:,2);
            
            errorL2  = sum((((uhq-uq).^2)*FE_SPACE.quad_weights').*MESH.jac');
            
            errorH1  = sum((((uhdxq-uxq).^2+(uhdyq-uyq).^2) *FE_SPACE.quad_weights').*MESH.jac');
            
            errorH1  = errorH1 + errorL2;
            
        end
        
    case 3
        
        x         = zeros(MESH.numElem,FE_SPACE.numQuadNodes);
        y         = x;
        z         = x;
        coord_ref = MESH.chi;
        
        for j = 1:4
            i = MESH.elements(j,:);
            vtemp = MESH.vertices(1,i);
            x = x + vtemp'*coord_ref(j,:);
            vtemp = MESH.vertices(2,i);
            y = y + vtemp'*coord_ref(j,:);
            vtemp = MESH.vertices(3,i);
            z = z + vtemp'*coord_ref(j,:);
        end
        
        uq      = DATA.uexact(x,y,z,t,DATA.param);
        
        uhe = uh(MESH.elements(1:FE_SPACE.numElemDof,:));
        uhq = uhe'*FE_SPACE.phi;
        
        if nargout == 1  % only L2 error is computed
            
            errorL2  = sum((((uhq-uq).^2)*FE_SPACE.quad_weights').*MESH.jac');
            
        elseif nargout == 2 % compute L2 and H1 errors
            
            uxq      = DATA.uxexact(x,y,z,t,DATA.param);
            uyq      = DATA.uyexact(x,y,z,t,DATA.param);
            uzq      = DATA.uzexact(x,y,z,t,DATA.param);
            
            dcdx     = MESH.invjac(:,1,1)';
            dedx     = MESH.invjac(:,1,2)';
            dtdx     = MESH.invjac(:,1,3)';
            dcdy     = MESH.invjac(:,2,1)';
            dedy     = MESH.invjac(:,2,2)';
            dtdy     = MESH.invjac(:,2,3)';
            dcdz     = MESH.invjac(:,3,1)';
            dedz     = MESH.invjac(:,3,2)';
            dtdz     = MESH.invjac(:,3,3)';
            
            dcdx     = dcdx(ones(FE_SPACE.numElemDof, 1), :);
            dedx     = dedx(ones(FE_SPACE.numElemDof, 1), :);
            dtdx     = dtdx(ones(FE_SPACE.numElemDof, 1), :);
            uhdxq    = (uhe.*dcdx)'*FE_SPACE.dphi_ref(:,:,1)+(uhe.*dedx)'*FE_SPACE.dphi_ref(:,:,2)+(uhe.*dtdx)'*FE_SPACE.dphi_ref(:,:,3);
            dcdy     = dcdy(ones(FE_SPACE.numElemDof, 1), :);
            dedy     = dedy(ones(FE_SPACE.numElemDof, 1), :);
            dtdy     = dtdy(ones(FE_SPACE.numElemDof, 1), :);
            uhdyq    = (uhe.*dcdy)'*FE_SPACE.dphi_ref(:,:,1)+(uhe.*dedy)'*FE_SPACE.dphi_ref(:,:,2)+(uhe.*dtdy)'*FE_SPACE.dphi_ref(:,:,3);
            dcdz     = dcdz(ones(FE_SPACE.numElemDof, 1), :);
            dedz     = dedz(ones(FE_SPACE.numElemDof, 1), :);
            dtdz     = dtdz(ones(FE_SPACE.numElemDof, 1), :);
            uhdzq    = (uhe.*dcdz)'*FE_SPACE.dphi_ref(:,:,1)+(uhe.*dedz)'*FE_SPACE.dphi_ref(:,:,2)+(uhe.*dtdz)'*FE_SPACE.dphi_ref(:,:,3);
            

            errorL2  = sum((((uhq-uq).^2)*FE_SPACE.quad_weights').*MESH.jac');
            
            errorH1  = sum((((uhdxq-uxq).^2+(uhdyq-uyq).^2+(uhdzq-uzq).^2) *FE_SPACE.quad_weights').*MESH.jac');
            
            errorH1  = errorH1 + errorL2;
            
        end
        
end

errorL2 = sqrt(errorL2);
errorH1 = sqrt(errorH1);


return