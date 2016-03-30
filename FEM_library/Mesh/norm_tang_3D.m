function [normal, normalf, FaceToElem_list, tang] = norm_tang_3D(boundaries, vertices, elements)
%NORM_TANG_3DP1 computes normal vectors on the boundary vertices of a P1 mesh
%
%   [NX,NY] = NORM_TANG_3DP1(BOUNDARIES,VERTICES,ELEMENTS) 
%   on the vertices of the boundaries of a P1 mesh compute the outward unit normal vectors 
%   NX, NY are arrrays of lenght NOV, where NOV Is the total number of mesh vertices.
%

%   Author: F. Negri (federico.negri@epfl.ch) 2013-2014
%   Copyright (C) Federico Negri, CMCS, EPFL
%


[n,nov]      = size(vertices);
[n,nside]    = size(boundaries);

[normalf, FaceToElem_list]  = faceNormal(vertices', boundaries', elements');
normalf                     = normalf';
normal                      = vertexNormal(vertices', boundaries', normalf')';

if nargout == 4
    t   = zeros(3,nov);
    % compute tangent vector
    for i = 1 : nov
        k   = find(normal(:,i),1);
        k_o = setdiff([1 2 3],k);
        t(k_o,i) = 1;
        t(k,i)   = (-t(:,i)'*normal(:,i))/normal(k,i);
        t(:,i)   = t(:,i)/norm(t(:,i));
    end
    
    tang = t;
end

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

%check correct orientation (outward)
faces    = faces';
elements = elements';
nodes = nodes';
elements_row = [elements(1,:) elements(2,:) elements(3,:) elements(4,:)];

id_e = zeros(size(faces,2),1);
parfor i = 1 : size(faces, 2)
    id_e(i)        = FaceToElement(faces(:,i), elements, elements_row);
    centroid_el    = mean(nodes(:,elements(:,id_e(i))),2) ;
    centroid_face  = mean(nodes(:,faces(:,i)),2) ;
    sign_normalf   = normals(i,:)*(centroid_face - centroid_el);
    if sign_normalf < 0
        normals(i,:) = - normals(i,:);
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

function normals = vertexNormal(vertices, faces, faceNormals)
%VERTEXNORMAL Compute normals to a mesh vertices
%
%   N = vertexNormal(V, F)
%   Computes vertex normals of the mesh given by vertices V and F. 
%   V is a vertex array with 3 columns, F is either a NF-by-3 or NF-by-4
%   index array, or a cell array with NF elements.
%
%   Example
%     % Draw the vertex normals of a sphere
%     s = [10 20 30 40];
%     [v f] = sphereMesh(s);
%     drawMesh(v, f);
%     view(3);axis equal; light; lighting gouraud;
%     normals = vertexNormal(v, f);
%     drawVector3d(v, normals);
%
%   See also
%     meshes3d, faceNormal, triangulateFaces
%
% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2011-12-19,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2011 INRA - Cepia Software Platform.


nv = size(vertices, 1);
nf = size(faces, 1);

% unit normals to the faces
if nargin < 3
    faceNormals = normalizeVector3d(faceNormal(vertices, faces));
end

% compute normal of each vertex: sum of normals to each face
normals = zeros(nv, 3);
if isnumeric(faces)
    for i = 1:nf
        face = faces(i, :);
        for j = 1:length(face)
            v = face(j);
            normals(v, :) = normals(v,:) + faceNormals(i,:);
        end
    end
else
    for i = 1:nf
        face = faces{i};
        for j = 1:length(face)
            v = face(j);
            normals(v, :) = normals(v,:) + faceNormals(i,:);
        end
    end
end

% normalize vertex normals to unit vectors
normals = normalizeVector3d(normals);
end
