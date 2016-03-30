function [nx,ny,tx,ty,normalf] = norm_tang_2D(boundaries,vertices,elements)
%NORM_TANG_2DP1 computes normal and tangential vectors on the boundary vertices of a P1 mesh
%
%   [NX,NY, TX,TY] = NORM_TANG_2DP1(BOUNDARIES,VERTICES,ELEMENTS) 
%   on the vertices of the boundaries of a P1 mesh compute the outward unit normal vectors 
%   and the corresponding tangential vectors. 
%   NX, NY, TX, TY are arrrays of lenght NOV, where NOV Is the total number of mesh vertices.
%
%   Author: F.Saleri, 12.04.04
%           F.Negri,  28.02.13
% -------------------------------
% Revision History


[n,nov]      = size(vertices);
[n,nside]    = size(boundaries);

nx           = zeros(nov,1);
ny           = zeros(nov,1);
tx           = nx;
ty           = ny;

elements_row = [elements(1,:) elements(2,:) elements(3,:)];
noe          = size(elements,2);

for iside = 1:nside
      i1 = boundaries(1,iside);
      i2 = boundaries(2,iside);
      l  = sqrt((vertices(1,i1)-vertices(1,i2))^2+(vertices(2,i1)-vertices(2,i2))^2);

      % coordinates of the vertices i1 and i2
      x1 = vertices(1,i1);
      x2 = vertices(1,i2);
      y1 = vertices(2,i1);
      y2 = vertices(2,i2);
      
      % normal of the edge
      n1 =   (y2-y1)/l;
      n2 = - (x2-x1)/l;
      
      % find the element to which the edge (i1,i2) belongs
      ie = find(elements_row == i1);
      for kk = 1:length(ie)  
            ie_tmp = ie(kk);
            if mod(ie_tmp,noe) == 0
                  ie_tmp = noe;
            else
                  ie_tmp = mod(ie_tmp,noe);
            end
            
            i3 = setdiff(elements(1:3,ie_tmp),[i1 i2]);
            if length(i3) == 1
                  break;
            end
      end
      
      % a rough check
      if length(i3) ~=1
            disp('error')
      end
      
      % coordinates of the third nodes of the elements to which 
      % the edge (i1,i2) belongs 
      x3  = vertices(1,i3);
      y3  = vertices(2,i3);

      % check if the normal is outward or not
      if (n1*(x3 - x2) + n2*(y3-y2)) > 0
            n1 = - n1;
            n2 = - n2;
      end
      
      normalf(:,iside) = [n1; n2];
      
      % normal and tangent vectors in the boundary vertices are computed as
      % the average of the corresponding vectors on the adjacent edges
      
      nx([i1 i2]) = nx([i1 i2]) + 0.5*n1;
      ny([i1 i2]) = ny([i1 i2]) + 0.5*n2;
      tx([i1 i2]) = tx([i1 i2]) + 0.5*n2;
      ty([i1 i2]) = ty([i1 i2]) + 0.5*n1;
      
end

return
