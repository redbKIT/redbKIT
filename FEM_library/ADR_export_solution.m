function ADR_export_solution(dim, u, vertices, elements, outputFileName, iter, variableName)
%ADR_EXPORT_SOLUTION export a (P1) finite element solution to vtk (binary) file
%
%   ADR_EXPORT_SOLUTION(DIM, U, VERTICES, ELEMENTS, OUTPUTFILENAME)
%   export the solution U to the vtk file OUTPUTFILENAME.vtk; DIM can be either
%   2 or 3, VERTICES is a DIM x numVertices matrix with the vertices
%   coordinates, ELEMENTS is the connectivity matrix. The default name of
%   the variable in the vtk file is 'u'
%
%   ADR_EXPORT_SOLUTION(DIM, U, VERTICES, ELEMENTS, OUTPUTFILENAME, ITER)
%   export the solution U to the vtk file OUTPUTFILENAME%04ITER
%
%   ADR_EXPORT_SOLUTION(DIM, U, VERTICES, ELEMENTS, OUTPUTFILENAME, ITER, VARIABLENAME)
%   the name of the variable in the vtk file is VARIABLENAME

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

titleData = 'Scalar_solution';
if nargin < 6
      iter = -1;
end

if nargin < 7
      variableName = 'u';
end

if dim == 2
      exportData=struct('iteration', {iter},...
            'vertices', {vertices'},...
            'elements', {elements(1:3,:)'},...
            'outputFile', {outputFileName},...
            'title', {titleData},...
            'variableName',{{variableName}},...
            'variableType',{{'SCALARS'}},...
            'variableData',{{u'}});
      
      exporter2dVTK(exportData);
      
elseif dim ==3
      exportData=struct('iteration', {iter},...
            'vertices', {vertices'},...
            'elements', {elements(1:4,:)'},...
            'outputFile', {outputFileName},...
            'title', {titleData},...
            'variableName',{{variableName}},...
            'variableType',{{'SCALARS'}},...
            'variableData',{{u'}});
      
      exporter3dVTK(exportData);
      
end
