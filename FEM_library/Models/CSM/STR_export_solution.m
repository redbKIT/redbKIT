function STR_export_solution(dim, u, vertices, elements, nov, outputFileName, iter, variableName)
%STR_EXPORT_SOLUTION export a (P1) vectorial finite element solution to vtk (binary) file
%
%   STR_EXPORT_SOLUTION(DIM, U, VERTICES, ELEMENTS, NOV, OUTPUTFILENAME)
%   export the solution U to the vtk file OUTPUTFILENAME.vtk; DIM can be either
%   2 or 3, VERTICES is a DIM x numVertices matrix with the vertices
%   coordinates, ELEMENTS is the connectivity matrix. The default name of
%   the variable in the vtk file is 'StructureDisplacement'
%
%   STR_EXPORT_SOLUTION(DIM, U, VERTICES, ELEMENTS, NOV, OUTPUTFILENAME, ITER)
%   export the solution U to the vtk file OUTPUTFILENAME%04ITER
%
%   STR_EXPORT_SOLUTION(DIM, U, VERTICES, ELEMENTS, NOV, OUTPUTFILENAME, ITER, VARIABLENAME)
%   the name of the variable in the vtk file is VARIABLENAME

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

titleData = 'STR_solution';
if nargin < 7
      iter = -1;
end

if nargin < 8
      variableName = 'StructureDisplacement';
end

novP1 = size(vertices,2);

if dim == 2
      
      exportData = struct('iteration', {iter},...
            'vertices', {vertices'},...
            'elements', {elements(1:3,:)'},...
            'outputFile', {outputFileName},...
            'title', {titleData},...
            'variableName',{{variableName}},...
            'variableType',{{'VECTORS'}},...
            'variableData',{{[full(u(1:novP1)) full(u(nov+1:novP1+nov))  0*full(u(1:novP1))]'}});
      
      exporter2dVTK(exportData);
      
elseif dim == 3
      
      exportData=struct('iteration', {iter},...
            'vertices', {vertices'},...
            'elements', {elements(1:4,:)'},...
            'outputFile', {outputFileName},...
            'title', {titleData},...
            'variableName',{{variableName}},...
            'variableType',{{'VECTORS'}},...
            'variableData',{{[full(u(1:novP1)) full(u(nov+1:novP1+nov))  full(u(2*nov+1:novP1+2*nov))]'}});
         
      exporter3dVTK(exportData);
end
