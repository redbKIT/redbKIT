function [elements,vertices,boundaries,rings]=P1toB1mesh2D(elements,vertices,boundaries,rings)
%P1TOB1MESH2D computes a B1 grid
%  [EB1,VB1,BP1]=P1TOP2MESH2D(EP1,VP1,BP1) computes a mesh for B1
%  finite elements starting from a P1 finite elements mesh.
%  EB1(1:3,:) is the connectivity matrix of the B1 mesh.
%  The number of the P1 vertices can be computed as max(max(EB1(1:3,:))).

%       F. Saleri 5-15-03.

[n,nov]=size(vertices);
[n,noe]=size(elements);

elements = [elements(1:3,:);[1:noe]+nov];
for ie = 1:noe
    xb = sum(vertices(1,elements(1:3,ie)))/3;
    yb = sum(vertices(2,elements(1:3,ie)))/3;
    vertices = [vertices, [xb;yb]];
end

    

return
