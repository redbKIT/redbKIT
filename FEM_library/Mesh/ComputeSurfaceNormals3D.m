function [normalf, FaceToElem_list] = ComputeSurfaceNormals3D(boundaries, vertices, elements)
%ComputeSurfaceNormals3D computes normal vectors on the boundary vertices 
%of a P1 TET mesh
%
%   [normalf] = ComputeSurfaceNormals3D(boundaries, vertices, elements)

%   This file is part of redbKIT.
%   Author: Federico Negri <federico.negri at epfl.ch>

if nargout == 1
    
    normalf  = faceNormal(vertices', boundaries', elements');
    
elseif nargout == 2
    
    [normalf, FaceToElem_list]  = faceNormal(vertices', boundaries', elements');
    
end

normalf  = normalf';

end



function [normals, id_e] = faceNormal(nodes, faces, elements)
%FACENORMAL Compute normal vector of faces in a 3D mesh
%
%   NORMALS = faceNormal(VERTICES, FACES)
%   VERTICES is a set of 3D points (as a N-by-3 array), and FACES is either
%   a N-by-3 index array or a cell array of indices. The function computes
%   the normal vector of each face.
%   The orientation of the normal is defined by the sign of cross product
%   between vectors joining vertices 1 to 2 and 1 to 3.
%
%
%   Example
%     [v e f] = createIcosahedron;
%     normals1 = faceNormal(v, f);
%     centros1 = faceCentroids(v, f);
%     figure; drawMesh(v, f); 
%     hold on; axis equal; view(3);
%     drawVector3d(centros1, normals1);
%
%     pts = rand(50, 3);
%     hull = minConvexHull(pts);
%     normals2 = faceNormal(pts, hull);
%
%   See also
%   meshes3d, drawMesh, convhull, convhulln, drawVector3d

% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2006-07-05
% Copyright 2006 INRA - CEPIA Nantes - MIAJ (Jouy-en-Josas).

% compute vector of first edges
v1 = nodes(faces(:,2),1:3) - nodes(faces(:,1),1:3);
v2 = nodes(faces(:,3),1:3) - nodes(faces(:,1),1:3);

% compute normals using cross product (nodes have same size)
normals = normalizeVector3d(cross(v1, v2, 2));

message = sprintf(['\n\n[\bWarning]\b: ComputeSurfaceNormals3D assumes that boundary elements have correct orientation,\n', ...
                  'so that surface normals are directed outward. No check on the sign is made.\n\n']);
fprintf(message);

%% check correct orientation (outward)
% Uncomment the following lines if your mesh does not have outward normals
%
% faces    = faces';
% elements = elements';
% nodes = nodes';
% elements_row = [elements(1,:) elements(2,:) elements(3,:) elements(4,:)];
% 
% id_e = zeros(size(faces,2),1);
% parfor i = 1 : size(faces, 2)
%     id_e(i)        = FaceToElement(faces(:,i), elements, elements_row);
%     centroid_el    = mean(nodes(:,elements(:,id_e(i))),2) ;
%     centroid_face  = mean(nodes(:,faces(:,i)),2) ;
%     sign_normalf   = normals(i,:)*(centroid_face - centroid_el);
%     if sign_normalf < 0
%         normals(i,:) = - normals(i,:);
%     end
% end

%% compute FaceToElem_list
if nargout > 1
    faces    = faces';
    elements = elements';
    elements_row = [elements(1,:) elements(2,:) elements(3,:) elements(4,:)];
    
    id_e = zeros(size(faces,2),1);
    parfor i = 1 : size(faces, 2)
        id_e(i)        = FaceToElement(faces(:,i), elements, elements_row);
    end
end
    
end

function id_e = FaceToElement(face, elements,   elements_row)
% find the element to which a face belongs
i1 = face(1);
i2 = face(2);
i3 = face(3);

noe = size(elements,2);

ie = find(elements_row == i1);
for kk = 1 : length(ie)

    ie_tmp = ie(kk);
    if mod(ie_tmp,noe) == 0
        ie_tmp = noe;%ie_tmp/noe;
    else
        ie_tmp = mod(ie_tmp,noe);
    end
    
    tmp = setdiff(elements(1:4,ie_tmp),[i1 i2 i3]);
    if length(tmp) == 1
        id_e  = ie_tmp;
        break;
    end
end

end

function vn = normalizeVector3d(v)
%NORMALIZEVECTOR3D Normalize a 3D vector to have norm equal to 1
%
%   V2 = normalizeVector3d(V);
%   Returns the normalization of vector V, such that ||V|| = 1. Vector V is
%   given as a row vector.
%
%   If V is a N-by-3 array, normalization is performed for each row of the
%   input array.
%
%   See also:
%   vectors3d, vectorNorm3d
%
%   ---------
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 29/11/2004.
%

% HISTORY
% 2005-11-30 correct a bug
% 2009-06-19 rename as normalizeVector3d
% 2010-11-16 use bsxfun (Thanks to Sven Holcombe)

vn   = bsxfun(@rdivide, v, sqrt(sum(v.^2, 2)));
end
