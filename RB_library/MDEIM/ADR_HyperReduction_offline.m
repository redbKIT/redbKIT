function [HyRED, U_matrix, U_rhs] = ADR_HyperReduction_offline(FOM, S_rhs, S_matrix, tolPOD)
%HYPERREDUCTION_OFFLINE given RHS and matrix snapshots, performs POD and
%(M)DEIM to generate an affine approximation of the system

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 


fig_folder = 'Figures/DEIM/';
[~,~,~] =  mkdir(fig_folder);

%% POD on RHS snapshots: U_rhs is a POD basis for the RHS
[U_rhs,    ~,    HyRED.Sigma_rhs]     = VPOD_basis_computation(S_rhs, [], tolPOD(1));

% generate figure
handle1 = figure;
subplot(1,2,1)
loglog(HyRED.Sigma_rhs./HyRED.Sigma_rhs(1),'-og');
hold on
cut_line_x = [1 : size(U_rhs,2)];
cut_line_y = HyRED.Sigma_rhs(size(U_rhs,2))./HyRED.Sigma_rhs(1) * ones(1,size(U_rhs,2));
loglog(cut_line_x, cut_line_y,'--k');
grid on
xlim([1  size(S_rhs,2)])
title('Rhs snapshots spectrum')

%% POD on matrix snapshots: U_matrix is a POD basis for the (vectorized) system matrix
[U_matrix,    ~,    HyRED.Sigma_matrix]     = MPOD_basis_computation(S_matrix, tolPOD(2));

% generate figure
subplot(1,2,2)
loglog(HyRED.Sigma_matrix./HyRED.Sigma_matrix(1),'-or');
hold on
cut_line_x = [1 : size(U_matrix,2)];
cut_line_y = HyRED.Sigma_matrix(size(U_matrix,2))./HyRED.Sigma_matrix(1) * ones(1,size(U_matrix,2));
loglog(cut_line_x, cut_line_y,'--k');
grid on
xlim([1  size(S_matrix,2)])
title('Matrix snapshots spectrum')

% save figure
saveas(handle1,strcat(fig_folder,'snapshots_spectrum'),'epsc');
saveas(handle1,strcat(fig_folder,'snapshots_spectrum'),'fig');

%% run (M)DEIM over matrices and RHS bases
[HyRED.IDEIM_rhs, ~, U_rhs]     = DEIM(U_rhs);
[HyRED.IDEIM_m ]                = MDEIM(U_matrix);

%% restrict POD bases to DEIM indices 
HyRED.U_matrix  = U_matrix(HyRED.IDEIM_m,:);
HyRED.U_rhs     = U_rhs(HyRED.IDEIM_rhs,:);
        
%% compute and store "reduced mesh" (corresponding to DEIM indices)
ndf       =  length(FOM.MESH.internal_dof); 

[ ~, node_to_element, node_to_boundary ] = compute_adjacency_elements(FOM.MESH.nodes, ...
    FOM.MESH.elements, FOM.MESH.dim, FOM.MESH.boundaries, FOM.FE_SPACE.fem); 

[IDEIM_mA_elem, ~, IDEIM_mA_bound ]     = ADR_DEIM_Index_to_Elements('matrix', HyRED.IDEIM_m,   ndf, node_to_element, ...
    node_to_boundary, FOM.MESH.internal_dof, FOM.MESH.numNodes);

[IDEIM_rhs_elem, ~, IDEIM_rhs_bound]    = ADR_DEIM_Index_to_Elements('rhs',    HyRED.IDEIM_rhs, ndf, node_to_element, ...
    node_to_boundary, FOM.MESH.internal_dof, FOM.MESH.numNodes);
  
IDEIM_all_elem       = unique([IDEIM_mA_elem  IDEIM_rhs_elem]);
IDEIM_all_bound      = unique([IDEIM_mA_bound  IDEIM_rhs_bound]);
IDEIM_all_nodes      = FOM.MESH.elements(:,IDEIM_all_elem);
IDEIM_all_nodes      = unique(IDEIM_all_nodes(:));

%% Save reduced mesh to vtk for visualization
ADR_export_solution(FOM.MESH.dim, ones(FOM.MESH.numVertices,1), FOM.MESH.vertices, FOM.MESH.elements(1:3,:), [fig_folder,'Reference_Mesh']);
ADR_export_solution(FOM.MESH.dim, ones(FOM.MESH.numVertices,1), FOM.MESH.vertices, FOM.MESH.elements(1:3,IDEIM_all_elem), [fig_folder,'Reduced_Mesh']);

%% Generate HyRed Structure

HyRED.Nelem           = length(IDEIM_all_elem);
HyRED.IDEIM_all_elem  = IDEIM_all_elem;
HyRED.IDEIM_all_bound = IDEIM_all_bound;
HyRED.IDEIM_all_nodes = IDEIM_all_nodes;

HyRED.MESH            = FOM.MESH;
HyRED.MESH.elements   = HyRED.MESH.elements(:,IDEIM_all_elem);
HyRED.MESH.numElem    = size(HyRED.MESH.elements,2);
HyRED.MESH.Neumann_side = intersect(HyRED.MESH.Neumann_side, IDEIM_all_bound);
HyRED.MESH.Robin_side   = intersect(HyRED.MESH.Robin_side,   IDEIM_all_bound);

HyRED.FE_SPACE         = FOM.FE_SPACE;
end