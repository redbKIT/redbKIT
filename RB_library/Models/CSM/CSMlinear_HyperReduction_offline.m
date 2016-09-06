function [HyRED, U_matrix, U_rhs] = CSMlinear_HyperReduction_offline(FOM, S_rhs, S_matrix, tolPOD)
%CSMlinear_HyperReduction_offline given RHS and matrix snapshots, performs POD and
%(M)DEIM to generate an affine approximation of the system

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 


fig_folder = 'Figures/DEIM/';
[~,~,~] =  mkdir(fig_folder);

%% POD on RHS snapshots: U_rhs is a POD basis for the RHS
[U_rhs,    ~,    HyRED.Sigma_rhs]     = VPOD_basis_computation(S_rhs, [], tolPOD(1), 1);

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
[U_matrix,    ~,    HyRED.Sigma_matrix]     = MPOD_basis_computation(S_matrix, tolPOD(2), 1);

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
RedMeshObject =  ReducedMesh( FOM.MESH, FOM.MESH.fem, 'CSM' );
RedMeshObject.AppendInternalDoFs_Vectorized( HyRED.IDEIM_m );
RedMeshObject.AppendInternalDoFs( HyRED.IDEIM_rhs );
RedMeshObject.Build( FOM.DATA );
RedMeshObject.ExportToVtk( fig_folder, 'CSMHyper');
Red_Mesh = RedMeshObject.M_Red_Mesh;

% for k = 1 : FOM.MESH.dim
%     Red_Mesh.Neumann_side{k} = intersect(Red_Mesh.Neumann_side{k}, RedMeshObject.M_ReducedBoundaries);
% end

%% Generate HyRed Structure
HyRED.MESH             = Red_Mesh;
HyRED.FE_SPACE         = FOM.FE_SPACE;

end