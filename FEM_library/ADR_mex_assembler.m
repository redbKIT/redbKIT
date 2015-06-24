function  [Arows, Acols, Acoef, Mcoef, Rrows, Rcoef] = ADR_mex_assembler(OPERATOR, TC_d, TC_t, elements, nln, mu, bx, by, si, f,...
    w,dcdx,dcdy,dedx,dedy,phi,dcsiphi,detaphi,detjac)
%ADR_MEX_ASSEMBLER core of the vectorized matrix assembly. 
% This function can be "mexified" (i.e. translated to C code, compiled 
% and linked via the mex interface ADR_mex_assembler_mex) using the
% command:
% coder -build ADR_assembly

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

OP = zeros(1,4);
switch OPERATOR
    
    case 'diffusion'
        OP(1) = 1;
    case 'transport'
        OP(2) = 1;
    case 'reaction'
        OP(3) = 1;
    case 'source'   
        OP(4) = 1;
    case 'all'
        OP    = ones(1,4);
end

[C_d, C_t] = set_TensorComponent(TC_d, TC_t);

%% Initialize vectors where to save (in "sparse format") entries and indices
% of the matrices and rhs
noe   = size(elements,2); % for parallel assembly

nln2  = nln*nln;
Arows = zeros(nln2*noe,1);
Acols = Arows;
Acoef = Arows;
Mcoef = Arows;
Rrows = zeros(nln*noe,1);
Rcoef = Rrows;

[rows,cols] = meshgrid(1:nln,1:nln);
rows        = rows(:);
cols        = cols(:);

iii         = 1:nln2;
ii          = 1:nln;

one         = ones(nln, 1);

%% Local mass matrix (computed only once)
% with quadrature nodes

MASS = (phi.*w(one,:))*phi';

%% Assembly: loop over the elements
for ie = 1 : noe
    muloc = mu(ie,:);
    bxloc = bx(ie,:);
    byloc = by(ie,:);
    muloc = muloc.*w;
    muloc = muloc(one,:);  % nln x nqn array
    bxloc = bxloc.*w;
    bxloc = bxloc(one,:);  % nln x nqn array
    byloc = byloc.*w;
    byloc = byloc(one,:);  % nln x nqn array
    siloc = si(ie,:).*w;
    siloc = siloc(one,:);  % nln x nqn array
    floc  = f(ie,:).*w;
    floc  = floc(1,:)';
    
    gradx = dcdx(ie)*dcsiphi + dedx(ie)*detaphi;
    grady = dcdy(ie)*dcsiphi + dedy(ie)*detaphi;
    
    trasp_loc    =  (C_t(1)*gradx.*bxloc + C_t(2)*grady.*byloc)*phi';
    
    react_loc    =  (phi.*siloc)*phi';
    
    diff_loc     =  C_d(1,1)*(gradx.*muloc)*gradx'   + C_d(1,2)*(gradx.*muloc)*grady' + ...
                    C_d(2,1)*(grady.*muloc)*gradx'   + C_d(2,2)*(grady.*muloc)*grady';
    
    aloc         =  OP(1)*diff_loc + OP(2)*trasp_loc + OP(3)*react_loc;
    
    Arows(iii)  =  elements(rows,ie);
    Acols(iii)  =  elements(cols,ie);
    Acoef(iii)  =  aloc(:)*detjac(ie);
    Mcoef(iii)  =  MASS(:)*detjac(ie);
    Rrows(ii)   =  elements(1:nln,ie);
    Rcoef(ii)   =  OP(4)*phi*floc*detjac(ie);
    iii         =  iii + nln2;
    ii          =  ii  + nln;
end

end


function [C_d, C_t] = set_TensorComponent(TC_d, TC_t)

if  TC_t(1)==10
    C_t      = [1 1];
else
    C_t          = zeros(1,2);
    C_t(TC_t(1)) = 1;
end


if (TC_d(1)==10 && TC_d(2)==10)
    C_d      = eye(2,2);
else
    C_d      = zeros(2,2);
    C_d(TC_d(1),TC_d(2)) = 1;
end

end