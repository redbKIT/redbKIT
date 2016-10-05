function [ FOM ] = build_affineFOM( elements, vertices, boundaries, fem, data_file, Advection_Field )
%   For reference, see Section 3.8, 6.6 and 7.2 of
%
%   Quarteroni, Manzoni, Negri - REDUCED BASIS METHODS FOR PARTIAL
%   DIFFERENTIAL EQUATIONS. AN INTRODUCTION. Springer, 2015
%
%   The high-fidelity (full-order) model is already provided in the FOM.mat
%   struct

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

%% For the moment HARDCODED inputs
FOM.Qa = 2;
FOM.Qf = 6;
FOM.P  = 4;
FOM.mu_min = [ 0  0  0  1 ];
FOM.mu_max = [ 12 12 12 600];

FOM.model = 'ADR';

%% Read problem parameters and BCs from data_file
DATA   = read_DataFile(data_file);
if nargin < 7
    DATA.param = [];
else
    DATA.param = param;
end
t      = [];

%% Set quad_order
quad_order                  = 5;  

%% Create and fill the MESH data structure
dim              = 3;
[ MESH ] = buildMESH( dim, elements, vertices, boundaries, fem, quad_order, DATA );

%% Create and fill the FE_SPACE data structure
[ FE_SPACE ] = buildFESpace( MESH, fem, 1, quad_order );

fprintf('\n **** PROBLEM''S SIZE INFO ****\n');
fprintf(' * Number of Vertices  = %d \n',MESH.numVertices);
fprintf(' * Number of Elements  = %d \n',MESH.numElem);
fprintf(' * Number of Nodes     = %d \n',MESH.numNodes);
fprintf(' * Number of Dofs      = %d \n',FE_SPACE.numDof);
fprintf('-------------------------------------------\n');


%% Fill FOM structures
FOM.MESH     = MESH;
FOM.FE_SPACE = FE_SPACE;
FOM.DATA     = DATA;

%% Assemble affine matrices and rhs vectors

fprintf('\n Assembling affine decomposition ...\n');


% A_1: diffusion
A_1       =  ADR_Assembler(MESH, DATA, FE_SPACE, 'diffusion', [], [], []);
FOM.Aq{1} =  A_1(FOM.MESH.internal_dof,FOM.MESH.internal_dof);

DATA.param  = [1 0 0];
[~, Fr]     = ADR_ApplyBC(A_1, [], FE_SPACE, MESH, DATA);
FOM.Fq{1}   = Fr;
DATA.param  = [0 1 0];
[~, Fr]     = ADR_ApplyBC(A_1, [], FE_SPACE, MESH, DATA);
FOM.Fq{2}   = Fr;
DATA.param  = [0 0 1];
[~, Fr]     = ADR_ApplyBC(A_1, [], FE_SPACE, MESH, DATA);
FOM.Fq{3}   = Fr;

% A_2: advection
DATA.transport{1} = Advection_Field(1:MESH.numNodes);
DATA.transport{2} = Advection_Field(1+MESH.numNodes:2*MESH.numNodes);
DATA.transport{3} = Advection_Field(1+2*MESH.numNodes:3*MESH.numNodes);

A_2       =  ADR_Assembler(MESH, DATA, FE_SPACE, 'transport', [], [], []);
FOM.Aq{2} =  A_2(FOM.MESH.internal_dof,FOM.MESH.internal_dof);

DATA.param  = [1 0 0];
[~, Fr]     = ADR_ApplyBC(A_2, [], FE_SPACE, MESH, DATA);
FOM.Fq{4}   = Fr;
DATA.param  = [0 1 0];
[~, Fr]     = ADR_ApplyBC(A_2, [], FE_SPACE, MESH, DATA);
FOM.Fq{5}   = Fr;
DATA.param  = [0 0 1];
[~, Fr]     = ADR_ApplyBC(A_2, [], FE_SPACE, MESH, DATA);
FOM.Fq{6}   = Fr;


FOM.u_D        =  @(x,mu)(0.*x(1,:).*x(2,:).*x(3,:) ...
               + mu(1) *  (x(2,:)>=-0.2) .*(x(2,:)<1) .* (x(1,:)<=1+0.1/2) .*(x(1,:)>=1-0.1/2) .*(x(3,:)>0) ... 
               + mu(2) *  (x(2,:)<=0.2)  .*(x(2,:)>-1) .* (x(1,:)<=2+0.1/2) .*(x(1,:)>=2-0.1/2) .*(x(3,:)>0)... 
               + mu(3) *  (x(2,:)>=-0.2) .*(x(2,:)<1) .* (x(1,:)<=3+0.1/2) .*(x(1,:)>=3-0.1/2).*(x(3,:)>0));


%% Compute Xnorm
DATA.diffusion =  @(x,y,z,t,param)(1+0.*x.*y);
X              =  ADR_Assembler(MESH, DATA, FE_SPACE, 'diffusion', [], [], []);
DATA.reaction  =  @(x,y,z,t,param)(1+0.*x.*y);
X              =  X + ADR_Assembler(MESH, DATA, FE_SPACE, 'reaction');

FOM.Xnorm =  X(FOM.MESH.internal_dof,FOM.MESH.internal_dof);
end
