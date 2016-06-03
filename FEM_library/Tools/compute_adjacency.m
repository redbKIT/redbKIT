function [ A ] = compute_adjacency(vertices, elements, dim, fem)
%COMPUTE_ADJACENCY compute adjacency matrix for 2d or 3d TRI/TET P1/P2 mesh

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

if nargin < 4 || isempty(fem)
    fem = 'P1';
end

noe = size(elements,2);
nov = size(vertices,2);

nln = select(fem, dim);

nln2 = nln^2;

[X,Y] = meshgrid(1:nln,1:nln);

rr = X(:)';
tt = Y(:)';
cc = ones(1,nln2);

row  = zeros(nln2*noe,1);
col  = row;
coef = row;

iii  = 1:nln2;

for i = 1 : noe
    
    p = elements(1:nln,i);
    
    row(iii)  = p(rr);
    col(iii)  = p(tt);
    coef(iii) = cc;
    
    iii = iii + nln2;
    
end

A = GlobalAssemble(row,col,coef,nov,nov);

end
