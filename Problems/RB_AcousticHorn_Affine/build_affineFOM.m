function [ FOM ] = build_affineFOM( elements, vertices, boundaries, fem, data_file )

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

%% For the moment HARDCODED inputs
FOM.Qa = 4;
FOM.Qf = 1;
FOM.P  = 1;
FOM.mu_min = [ 10    ];
FOM.mu_max = [ 1800  ];
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
    FOM.Aq{i} = sparse(length(FOM.MESH.internal_dof), length(FOM.MESH.internal_dof));
end

for i = 1 : FOM.Qf
    FOM.Fq{i} = sparse(length(FOM.MESH.internal_dof), 1);
end

% A_1: diffusion
A_1       =  ADR_Assembler(MESH, DATA, FE_SPACE, 'diffusion');
FOM.Aq{1} =  A_1(FOM.MESH.internal_dof,FOM.MESH.internal_dof);

% A_2: reaction
A_2       =  ADR_Assembler(MESH, DATA, FE_SPACE, 'reaction');
FOM.Aq{2} =  A_2(FOM.MESH.internal_dof,FOM.MESH.internal_dof);

% A_3 : mass matrix on portions \Gamma_o (flag 3) and \Gamma_i (flag 2) of Robin boundaries
flag_Gamma_o    = 3;
flag_Gamma_i    = 2;
GAMMA_side_o    = find(FOM.MESH.boundaries(5,:)==flag_Gamma_o);
GAMMA_side_i    = find(FOM.MESH.boundaries(5,:)==flag_Gamma_i);
GAMMA_side      = unique([GAMMA_side_o GAMMA_side_i]);

GAMMA_vert_oi = FOM.MESH.boundaries(1:FOM.MESH.numBoundaryDof, GAMMA_side);
GAMMA_vert_oi = unique(GAMMA_vert_oi(:)); 

[~, A_3] = AssembleMass1D( fem, FOM.MESH.boundaries, FOM.MESH.vertices, FOM.MESH.numNodes, GAMMA_side, GAMMA_vert_oi);

FOM.Aq{3} =  A_3(FOM.MESH.internal_dof,FOM.MESH.internal_dof);

% A_4 : mass matrix on portion \Gamma_o (flag 3) of Robin boundaries
GAMMA_vert_o = FOM.MESH.boundaries(1:FOM.MESH.numBoundaryDof, GAMMA_side_o);
GAMMA_vert_o = unique(GAMMA_vert_o(:));

[~, A_4]    = AssembleMass1D( fem, FOM.MESH.boundaries, FOM.MESH.vertices, FOM.MESH.numNodes, GAMMA_side_o, GAMMA_vert_o);
FOM.Aq{4}   =  A_4(FOM.MESH.internal_dof,FOM.MESH.internal_dof);

% F_1: diffusion
GAMMA_vert_i = FOM.MESH.boundaries(1:FOM.MESH.numBoundaryDof, GAMMA_side_i);
GAMMA_vert_i = unique(GAMMA_vert_i);
[~, M_i]     = AssembleMass1D( fem, FOM.MESH.boundaries, FOM.MESH.vertices, FOM.MESH.numNodes, GAMMA_side_i, GAMMA_vert_i);
u_i = zeros(length(FOM.MESH.internal_dof), 1);
u_i(GAMMA_vert_i) = 1;
FOM.Fq{1}      =  M_i * u_i;

FOM.u_D        =  @(x,mu)(FOM.DATA.bcDir(x(1,:),x(2,:),[],mu));


%% Compute Xnorm
DATA.diffusion =  @(x,y,t,param)(1+0.*x.*y);
X              =  ADR_Assembler(MESH, DATA, FE_SPACE, 'diffusion', [], [], []);
DATA.reaction  =  @(x,y,t,param)(1+0.*x.*y);
M              =  ADR_Assembler(MESH, DATA, FE_SPACE, 'reaction');

FOM.Xnorm =  X(FOM.MESH.internal_dof,FOM.MESH.internal_dof) + M(FOM.MESH.internal_dof,FOM.MESH.internal_dof);
end
