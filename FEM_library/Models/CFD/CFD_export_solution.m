function CFD_export_solution(dim, u, p, vertices, elements, numDofsVel, outputFileName, iter, variableName)
%CFD_EXPORT_SOLUTION export a (P1) vectorial finite element solution to vtk (binary) file
%
%   CFD_EXPORT_SOLUTION(DIM, U, P, VERTICES, ELEMENTS, NUMDOFSVEL, OUTPUTFILENAME)
%   export the solution U, P to the vtk file OUTPUTFILENAME.vtk; DIM can be either
%   2 or 3, VERTICES is a DIM x numVertices matrix with the vertices
%   coordinates, ELEMENTS is the connectivity matrix. The default name of
%   the variable in the vtk file is {'Pressure','Velocity'}
%
%   CFD_EXPORT_SOLUTION(DIM, U, P, VERTICES, ELEMENTS, NUMDOFSVEL, OUTPUTFILENAME, ITER)
%   export the solution U to the vtk file OUTPUTFILENAME%04ITER
%
%   CFD_EXPORT_SOLUTION(DIM, U, P, VERTICES, ELEMENTS, NUMDOFSVEL, OUTPUTFILENAME, ITER, VARIABLENAME)
%   the name of the variable in the vtk file is VARIABLENAME

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

titleData = 'STR_solution';
if nargin < 8
      iter = -1;
end

if nargin < 9
      variableName = {'Pressure', 'Velocity'};
end

novP1 = size(vertices,2);

if dim == 2
      
      exportData = struct('iteration', {iter},...
            'vertices', {vertices'},...
            'elements', {elements(1:3,:)'},...
            'outputFile', {outputFileName},...
            'title', {titleData},...
            'variableName',{variableName},...
            'variableType',{{'SCALARS', 'VECTORS'}},...
            'variableData',{{full(p') [full(u(1:novP1)) full(u(numDofsVel+1:novP1+numDofsVel))  0*full(u(1:novP1))]'}});
      
      exporter2dVTK(exportData);
      
elseif dim == 3
      
      exportData=struct('iteration', {iter},...
            'vertices', {vertices'},...
            'elements', {elements(1:4,:)'},...
            'outputFile', {outputFileName},...
            'title', {titleData},...
            'variableName',{variableName},...
            'variableType',{{'SCALARS', 'VECTORS'}},...
            'variableData',{{full(p') [full(u(1:novP1)) full(u(numDofsVel+1:novP1+numDofsVel))  full(u(2*numDofsVel+1:novP1+2*numDofsVel))]'}});
         
      exporter3dVTK(exportData);
end
