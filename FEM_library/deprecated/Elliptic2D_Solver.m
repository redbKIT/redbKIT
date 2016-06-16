function [u, FE_SPACE, MESH, DATA, errorL2, errorH1] = Elliptic2D_Solver(varargin)
%ELLIPTIC2D_SOLVER 2D diffusion-transport-reaction finite element solver
%
%   [U, FE_SPACE, MESH, DATA, ERRORL2, ERRORH1] = ...
%    ELLIPTIC2D_SOLVER(ELEMENTS, VERTICES, BOUNDARIES, FEM, DATA_FILE, PARAM)
%
%   Inputs:
%     ELEMENTS, VERTICES, BOUNDARIES: mesh information
%     FEM: string 'P1' or 'P2'
%     DATA_FILE: name of the file defining the problem data and
%          boundary conditions.
%     PARAM: vector of parameters possibly used in the data_file; 
%         if not provided, the PARAM vector is set to the empty vector.
%
%   Outputs:
%     U: problem solution
%     ERRORL2: L2-error between the numerical solution and the exact one 
%        (provided by the user in the data_file)
%     ERRORH1: H1-error between the numerical solution and the exact one 
%        (provided by the user in the data_file)
%     FE_SPACE: struct containing Finite Element Space information
%     MESH: struct containing mesh information
%     DATA: struct containing problem data

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 


warning('Elliptic2D_Solver is deprecated. Use Elliptic_Solver instead.')

[u, FE_SPACE, MESH, DATA, errorL2, errorH1] = Elliptic_Solver(2, varargin{:});


end
