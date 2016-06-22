function [elements, vertices, boundaries, rings] = P1toP2mesh2D(elements,vertices,boundaries, rings)
%P1TOP2MESH1D builds a P2 mesh in 2D.
%       F. Saleri 9-20-01.
%       F. Negri 2016, Add mesh Graph to speedup computations (still inefficient)

[~,nov]     =  size(vertices);
[~,noe]     =  size(elements);
nside       =  nov;

elements = [elements(1:3,:); zeros(3,noe); elements(4,:)];

[ a ] = compute_adjacency(vertices, elements, 2, 'P1');

[ii,jj,vv] = find(a);
a = sparse(ii, jj, vv*0 - 1, nov, nov);

for ie = 1:noe
    i = elements(1,ie);
    j = elements(2,ie);
    k = elements(3,ie);
    
    l1 = a(i,j);
    if l1 == -1 
       nside = nside + 1;
       a(i,j) = nside;
       a(j,i) = nside;
       elements(4,ie) = nside;
       vertices(1,nside) = (vertices(1,i)+vertices(1,j))*0.5;
       vertices(2,nside) = (vertices(2,i)+vertices(2,j))*0.5;
   else 
       elements(4,ie) = l1;
   end
   
   l2 = a(j,k);
   if l2 == -1 
       nside = nside + 1;
       a(j,k) = nside;
       a(k,j) = nside;
       elements(5,ie) = nside;
       vertices(1,nside) = (vertices(1,j)+vertices(1,k))*0.5;
       vertices(2,nside) = (vertices(2,j)+vertices(2,k))*0.5;
   else 
       elements(5,ie) = l2;
   end
   
   l3 = a(k,i);
   if l3 == -1 
       nside = nside + 1;
       a(k,i) = nside;
       a(i,k) = nside;
       elements(6,ie) = nside;
       vertices(1,nside) = (vertices(1,k)+vertices(1,i))*0.5;
       vertices(2,nside) = (vertices(2,k)+vertices(2,i))*0.5;
   else
       elements(6,ie) = l3;
   end
end

[n,nside] = size(boundaries);
for i = 1 : nside
    boundaries(3,i) = a(boundaries(1,i),boundaries(2,i));
end


return