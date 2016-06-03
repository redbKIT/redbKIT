function [ FOM ] = build_affineFOM( elements, vertices, boundaries, fem, data_file )
%   mu1: x-coordinate of the center of the Gaussian
%   mu2: y-coordinate of the center of the Gaussian
%   mu3: angle of the trasport field

%   For reference, see Section 10.5 of
%
%   Quarteroni, Manzoni, Negri - REDUCED BASIS METHODS FOR PARTIAL
%   DIFFERENTIAL EQUATIONS. AN INTRODUCTION. Springer, 2015

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 


%% For the moment HARDCODED inputs
FOM.Qa = 3;
FOM.P  = 3;
FOM.mu_min =  [     0.2     0.15   0  ];
FOM.mu_max =  [     0.8     0.35   360  ];
FOM.model  = 'ADR';

%% Read problem parameters and BCs from data_file
DATA       = read_DataFile(data_file);
DATA.param = [];
t          = [];

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

%% Assemble affine matrices

% initialize
for i = 1 : FOM.Qa
    FOM.Aq{i} = sparse(FOM.FE_SPACE.numDof, FOM.FE_SPACE.numDof);
end


% A_1: diffusion + reaction
A_1d       =  ADR_Assembler(MESH, DATA, FE_SPACE, 'diffusion');
A_1r       =  ADR_Assembler(MESH, DATA, FE_SPACE, 'reaction');
FOM.Aq{1}  =  A_1d(FOM.MESH.internal_dof, FOM.MESH.internal_dof) ...
            + A_1r(FOM.MESH.internal_dof, FOM.MESH.internal_dof);

% A_2: transport
A_2       =  ADR_Assembler(MESH, DATA, FE_SPACE, 'transport', [], 1);
FOM.Aq{2} =  A_2(FOM.MESH.internal_dof,FOM.MESH.internal_dof);

% A_3: transport
A_3       =  ADR_Assembler(MESH, DATA, FE_SPACE, 'transport', [], 2);
FOM.Aq{3} =  A_3(FOM.MESH.internal_dof,FOM.MESH.internal_dof);


%% EIM on RHS
reset_randomSeed;
Ns         =  1000;
EIM_mu_min =  [ 0.2     0.15     ];
EIM_mu_max =  [ 0.8     0.35     ];
mu_cube    =  lhsdesign(Ns, 2); % generate normalized design
mu_EIM     =  bsxfun(@plus,EIM_mu_min,bsxfun(@times,mu_cube,(EIM_mu_max-EIM_mu_min)));

x = zeros(MESH.numElem,FE_SPACE.numQuadNodes); y = x;
coord_ref = MESH.chi;

for j = 1 : 3
      i = MESH.elements(j,:);
      vtemp = MESH.vertices(1,i);
      x = x + vtemp'*coord_ref(j,:);
      vtemp = MESH.vertices(2,i);
      y = y + vtemp'*coord_ref(j,:);
end
quad_nodes = [x(:) y(:)];
 
tol_EIM            = 1e-3;
Max_Iter           = 100;
[IEIM, PHI_EIM]    = EmpiricalInterpolation(@(x,mu)nonaffine_source(x,mu), quad_nodes, mu_EIM, tol_EIM, Max_Iter);
FOM.Qf             = length(IEIM);

for i = 1 : FOM.Qf
    fprintf('Assemble source term #%d \n', i)
    DATA.force     =  reshape(PHI_EIM(:,i), MESH.numElem, FE_SPACE.numQuadNodes);
    [~, F_1]       =  ADR_Assembler(MESH, DATA, FE_SPACE, 'source');
    [~, F_1]       =  ADR_ApplyBC([], F_1, FE_SPACE, MESH, DATA);
    FOM.Fq{i}      =  F_1;
end

PHI_EIM  = PHI_EIM(IEIM, :);
EIM_quad = quad_nodes(IEIM, :);

[~,~,~] = mkdir('EIM');
save('EIM/IEIM.mat', 'IEIM');
save('EIM/EIM_quad.mat', 'EIM_quad');
save('EIM/PHI_EIM.mat', 'PHI_EIM');
addpath([pwd '/EIM'])

FOM.u_D        =  @(x, mu)(FOM.DATA.bcDir(x(1,:),x(2,:),[],mu));

%% Compute Xnorm
DATA.diffusion =  @(x,y,t,param)(1+0.*x.*y);
X              =  ADR_Assembler(MESH, DATA, FE_SPACE, 'diffusion', [], [], []);
DATA.reaction  =  @(x,y,t,param)(1+0.*x.*y);
M              =  ADR_Assembler(MESH, DATA, FE_SPACE, 'reaction');

FOM.Xnorm =  X(FOM.MESH.internal_dof,FOM.MESH.internal_dof) +...
             M(FOM.MESH.internal_dof,FOM.MESH.internal_dof);
end
