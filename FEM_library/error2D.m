function [errorL2,errorH1] = error2D(uh, MESH, DATA, FE_SPACE, t)
%ERROR2D compute L2, H1 and RHO error for 2D fem approximation
%
%   [ERRORL2] = ERROR2D(UH, MESH, DATA, FE_SPACE) computes
%   the L^2 error for FEM finite elements solution UH. The exact solution 
%   must be precised in the DATA.uexact field (see datafile template). 
%
%   [ERRORL2,ERRORH1] = ERROR2D(UH, MESH, DATA, FE_SPACE) computes
%   the L^2 and the H^1 error. In this case the struct DATA must provide
%   also the field uxexact and uyexact for the first derivatives 
%   of the exact solution.
%
%   [ERRORL2,ERRORH1] = ERROR1D(UH, MESH, DATA, FE_SPACE, T) computes
%   the errors at time T for time-dependent problems.

%   Auhtors: F. Saleri 09-01-05. F. Negri 18.11.2014

if nargin == 4
    t = [];
end

errorL2  = 0;
errorH1  = 0;

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
    uhdxq    = (uhe.*dcdx)'*FE_SPACE.dcsiphi+(uhe.*dedx)'*FE_SPACE.detaphi;
    dcdy     = dcdy(ones(FE_SPACE.numElemDof, 1), :);     
    dedy     = dedy(ones(FE_SPACE.numElemDof, 1), :);     
    uhdyq    = (uhe.*dcdy)'*FE_SPACE.dcsiphi+(uhe.*dedy)'*FE_SPACE.detaphi;
    
    errorL2  = sum((((uhq-uq).^2)*FE_SPACE.quad_weights').*MESH.jac');
    
    errorH1  = sum((((uhdxq-uxq).^2+(uhdyq-uyq).^2) *FE_SPACE.quad_weights').*MESH.jac');
     
    errorH1  = errorH1 + errorL2;

end

errorL2 = sqrt(errorL2);
errorH1 = sqrt(errorH1);


return