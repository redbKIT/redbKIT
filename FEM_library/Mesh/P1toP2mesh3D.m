function [elements,vertices,boundaries,rings]=P1toP2mesh3D(elements,vertices,boundaries,rings)
%P1TOP2MESH3D computes a P2 grid in 3D

%   F. Saleri 9-20-01.
%   F. Negri 2016, Add mesh Graph to speedup computations (still inefficient)

[~,nov]     =  size(vertices);
[~,noe]     =  size(elements);
nside       =  nov;

elements = [elements(1:4,:); zeros(6,noe); elements(5,:)];
 
[ a ] = compute_adjacency(vertices, elements, 3, 'P1');

[ii,jj,vv] = find(a);
a = sparse(ii, jj, vv*0 - 1, nov, nov);

for ie = 1:noe
      i = elements(1,ie);%0
      j = elements(2,ie);%1
      k = elements(3,ie);%2
      h = elements(4,ie);%3
      
      % 4, mid point of 0-1
      l1 = a(i,j);
      if l1 == -1
            nside = nside + 1;
            a(i,j) = nside;
            a(j,i) = nside;
            elements(5,ie) = nside;
            vertices(1:3,nside) = (vertices(1:3,i)+vertices(1:3,j))*0.5;
      else
            elements(5,ie) = l1;
      end
      
      % 5, mid point of 1-2
      l2 = a(j,k);
      if l2 == -1
            nside = nside + 1;
            a(j,k) = nside;
            a(k,j) = nside;
            elements(6,ie) = nside;
            vertices(1:3,nside) = (vertices(1:3,j)+vertices(1:3,k))*0.5;
      else
            elements(6,ie) = l2;
      end
      
      % 6, mid point of 2-0
      l3 = a(k,i);
      if l3 == -1
            nside = nside + 1;
            a(k,i) = nside;
            a(i,k) = nside;
            elements(7,ie) = nside;
            vertices(1:3,nside) = (vertices(1:3,k)+vertices(1:3,i))*0.5;
      else
            elements(7,ie) = l3;
      end
      
      % 7, mid point of 0-3
      l4 = a(i,h);
      if l4 == -1
            nside = nside + 1;
            a(h,i) = nside;
            a(i,h) = nside;
            elements(8,ie) = nside;
            vertices(1:3,nside) = (vertices(1:3,h)+vertices(1:3,i))*0.5;
      else
            elements(8,ie) = l4;
      end
      
      % 8, mid point of 1-3
      l5 = a(j,h);
      if l5 == -1
            nside = nside + 1;
            a(h,j) = nside;
            a(j,h) = nside;
            elements(9,ie) = nside;
            vertices(1:3,nside) = (vertices(1:3,h)+vertices(1:3,j))*0.5;
      else
            elements(9,ie) = l5;
      end
      
      % 9, mid point of 2-3
      l6 = a(k,h);
      if l6 == -1
            nside = nside + 1;
            a(h,k) = nside;
            a(k,h) = nside;
            elements(10,ie) = nside;
            vertices(1:3,nside) = (vertices(1:3,k)+vertices(1:3,h))*0.5;
      else
            elements(10,ie) = l6;
      end
      
     
end

[n,nside]=size(boundaries);
for i = 1 : nside
      boundaries(4,i) = a(boundaries(1,i),boundaries(2,i));
      boundaries(5,i) = a(boundaries(2,i),boundaries(3,i));
      boundaries(6,i) = a(boundaries(3,i),boundaries(1,i));
end

if nargin == 4
    
    [~,nrings]=size(rings);
    for i = 1 : nrings
        rings(3,i) = a(rings(1,i),rings(2,i));
    end
    
else
    rings = [];
end

return
