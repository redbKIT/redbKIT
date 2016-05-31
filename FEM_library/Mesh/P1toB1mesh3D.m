function [elements,vertices,boundaries,rings]=P1toB1mesh3D(elements,vertices,boundaries,rings)
%P1TOB1MESH3D build a B1 mesh in 3D

%   Author: F. Negri (federico.negri@epfl.ch) 2013-2014
%   Copyright (C) Federico Negri, CMCS, EPFL

[n,nov]  =  size(vertices);
[n,noe]  =  size(elements);

elements = [elements(1:4,:); [1:noe]+nov; elements(5:end,:)];

for ie = 1:noe
    
    xb = sum(vertices(1,elements(1:4,ie)))/4;
    yb = sum(vertices(2,elements(1:4,ie)))/4;
    zb = sum(vertices(3,elements(1:4,ie)))/4;
    vertices = [vertices, [xb;yb;zb] ];
    
end



return
