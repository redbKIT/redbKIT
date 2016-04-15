function [ FOM ] = build_basicFOM(elements, vertices, boundaries, fem, data_file)
%BUILD_AFFINEFOM returns a FOM (for full-order model) struct containing all
%the structures required for the reduction.

%   Author: F. Negri (federico.negri@epfl.ch) 2015
%   Copyright (C) Federico Negri, CMCS, EPFL

%% Parameters
FOM.P  = 5;
FOM.mu_min = [ 50    -0.03     -0.03     -0.03    -0.03];
FOM.mu_max = [ 1000  0.03       0.03      0.03     0.03];
FOM.mu_bar = [ 500    0          0         0        0]; 

%% Read problem parameters and BCs from data_file
DATA       = read_DataFile(data_file);
DATA.param = zeros(1,5);
DATA.shape_param = [2:5];
t          = [];
%% Set quad_order
quad_order                  = 4;  

%% Fill MESH data structure
dim              = 2;
[ MESH ] = buildMESH( dim, elements, vertices, boundaries, fem, quad_order, DATA );

%% Create and fill the FE_SPACE data structure
[ FE_SPACE ] = buildFESpace( MESH, fem, 1, quad_order );

%% Fill FOM structures
FOM.MESH     = MESH;
FOM.FE_SPACE = FE_SPACE;
FOM.DATA     = DATA;

%% Compute Xnorm
DATA.diffusion =  @(x,y,t,param)(1+0.*x.*y);
X              =  ADR_Assembler(MESH, DATA, FE_SPACE, 'diffusion', [], [], []);
DATA.reaction  =  @(x,y,t,param)(1+0.*x.*y);
M              =  ADR_Assembler(MESH, DATA, FE_SPACE, 'reaction');

FOM.Xnorm =  X(FOM.MESH.internal_dof,FOM.MESH.internal_dof) +...
             M(FOM.MESH.internal_dof,FOM.MESH.internal_dof);

 
end
