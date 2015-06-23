%READ_MSH reads a mesh in  MSH (ASCII) format version 2
%
%   MESH = READ_MSH(FILENAME) given a string containing the name of a msh
%   file,  returns a MESH struct containing the following fields:
%
%   MESH.elm_types: a 1x31 struct array with fields 
%                   NumNodesPerElem, Description, HowMany  
%
%   MESH.NODES: a 3xNumNodes matrix containing the nodes coordinates
%
%   MESH.ELEMENTS: a 31x1 cell array. The i-th cell contains a connectivity
%      matrix of size NumNodesPerElem x (HowMany+1) of elements of type i; 
%      the first NumNodesPerElem rows contain the IDs of the nodes belonging 
%      to the elements; the last row contains the the physical entity to 
%      which the elements belong 
%
%   For a description of the MSH ASCII File Format see 
%   http://www.manpagez.com/info/gmsh/gmsh-2.2.6/gmsh_63.php
%   or
%   http://geuz.org/gmsh/doc/texinfo/gmsh.html#MSH-ASCII-file-format
