function [ FOM ] = build_affineFOM( elements, vertices, boundaries, fem, data_file )
%   For reference, see Section 8.3 and 9.1 of
%
%   Quarteroni, Manzoni, Negri - REDUCED BASIS METHODS FOR PARTIAL
%   DIFFERENTIAL EQUATIONS. AN INTRODUCTION. Springer, 2015

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

%% For the moment HARDCODED inputs
FOM.Qa = 6;
FOM.Qf = 1;
FOM.P  = 3;
FOM.mu_min = [ 0.0      1      2     ];
FOM.mu_max = [ 0.6      20     50    ];

FOM.model = 'ADR';

%% Read problem parameters and BCs from data_file
DATA   = read_DataFile(data_file);
if nargin < 6
    DATA.param = [];
else
    DATA.param = param;
end
t      = [];

%% Set quad_order
quad_order                  = 4;  

%% Create and fill the MESH data structure
dim              = 2;
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

% initialize
for i = 1 : FOM.Qa
    FOM.Aq{i} = sparse(FOM.FE_SPACE.numDof, FOM.FE_SPACE.numDof);
end

for i = 1 : FOM.Qf
    FOM.Fq{i} = sparse(FOM.FE_SPACE.numDof, 1);
end

% A_1: diffusion
A_1       =  ADR_Assembler(MESH, DATA, FE_SPACE, 'diffusion', [1 1], [], 1);
FOM.Aq{1} =  A_1(FOM.MESH.internal_dof,FOM.MESH.internal_dof);

% A_2: diffusion
A_2       =  ADR_Assembler(MESH, DATA, FE_SPACE, 'diffusion', [2 2], [], 1);
FOM.Aq{2} =  A_2(FOM.MESH.internal_dof,FOM.MESH.internal_dof);


% A_3: transport
A_3       =  ADR_Assembler(MESH, DATA, FE_SPACE, 'transport', [], 2, 1);
FOM.Aq{3} =  A_3(FOM.MESH.internal_dof,FOM.MESH.internal_dof);

% A_4: diffusion
A_4       =  ADR_Assembler(MESH, DATA, FE_SPACE, 'diffusion', [], [], 2);
FOM.Aq{4} =  A_4(FOM.MESH.internal_dof,FOM.MESH.internal_dof);

% A_5: diffusion
A_5       =  ADR_Assembler(MESH, DATA, FE_SPACE, 'diffusion', [], [], 3);
FOM.Aq{5} =  A_5(FOM.MESH.internal_dof,FOM.MESH.internal_dof);

% A_6: diffusion
A_6       =  ADR_Assembler(MESH, DATA, FE_SPACE, 'diffusion', [], [], 4);
FOM.Aq{6} =  A_6(FOM.MESH.internal_dof,FOM.MESH.internal_dof);

% F_1: diffusion
[~, F_1]       =  ADR_Assembler(MESH, DATA, FE_SPACE, 'source', [], [], 2);
[~, F_1, u_D]  =  ADR_ApplyBC([], F_1, FE_SPACE, MESH, DATA);
FOM.Fq{1}      =  F_1;

FOM.u_D        =  @(x,mu)(FOM.DATA.bcDir(x(1,:),x(2,:),[],mu));


%% Compute Xnorm
DATA.diffusion =  @(x,y,t,param)(1+0.*x.*y);
X              =  ADR_Assembler(MESH, DATA, FE_SPACE, 'diffusion', [], [], []);
DATA.reaction  =  @(x,y,t,param)(1+0.*x.*y);
M              =  ADR_Assembler(MESH, DATA, FE_SPACE, 'reaction');

FOM.Xnorm =  X(FOM.MESH.internal_dof,FOM.MESH.internal_dof);
end
