function [] = test()

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>   

clc
clear all

[~,~,~] = mkdir('SnapshotsTest');
[~,~,~] = mkdir('FiguresTest');

delete('SnapshotsTest/SystemSnapshots.h5')
delete('SnapshotsTest/DisplacementSnapshots.h5')

%% ========================================================================
% DATA
dim      =  3;
fem      =  'P1';

% Parameters range
P = 3;% Yuoung modulus, Poisson Coefficient, External load
mu_min = [ 6*10^4   0.3  1000];
mu_max = [ 7*10^4   0.4  2000];
mu_bar = [ 6.5*10^4 0.35 1500]; 
    
% Training Parameters
mu_train_Dimension   = 10; % number of samples
mu_cube              = lhsdesign(mu_train_Dimension,P); % generate normalized design
Training_Parameters  = bsxfun(@plus,mu_min,bsxfun(@times,mu_cube,(mu_max-mu_min)));

OfflineTraining.Solution.h5_filename       = 'SnapshotsTest/DisplacementSnapshots.h5';
OfflineTraining.Solution.h5_section        = 'Displacement';

OfflineTraining.System.h5_filename                      = 'SnapshotsTest/SystemSnapshots.h5';
OfflineTraining.System.InternalForces.h5_section        = 'F_int';
OfflineTraining.System.ExternalForces.h5_section        = 'F_ext';

tol_POD_U    = 1e-4;
tol_POD_Fext = 1e-4;
tol_POD_Fint = 1e-4;
%% ========================================================================
% Solve High-Fidelity Models and Collect solution Snapshots

% load P1 mesh
[vertices, boundaries, elements] = msh_to_Mmesh('Cube', dim);

for i = 1 : mu_train_Dimension

    CSM_Solver(dim, elements, vertices, boundaries, fem, 'datafileR', ...
               Training_Parameters(i,:), ['FiguresTest/TrainingTierI_', num2str(i)], [], OfflineTraining.Solution);

end

%% ========================================================================
% Generate POD basis for the solution
DispSnap_HDF5  = HDF5_DenseMultiCVector(OfflineTraining.Solution.h5_filename, OfflineTraining.Solution.h5_section);
S_u            = DispSnap_HDF5.readValues();
[V,    ~,    Sigma_sol]  = VPOD_basis_computation(S_u, [], tol_POD_U, 1);

figure
loglog(Sigma_sol./Sigma_sol(1),'-or');
hold on
cut_line_x = [1 : size(V,2)];
cut_line_y = Sigma_sol(size(V,2))./Sigma_sol(1) * ones(1,size(V,2));
loglog(cut_line_x, cut_line_y,'--k');
grid on
xlim([1  size(S_u,2)])
title('Solution snapshots spectrum')

clear S_u;

ROM.V = V;

%% ========================================================================
% Solve POD-Galerkin ROM and Collect System Snapshots

for i = 1 : mu_train_Dimension

    CSM_POD_Solver(dim, elements, vertices, boundaries, fem, 'datafileR', ...
               Training_Parameters(i,:), ['FiguresTest/TrainingTierII_', num2str(i)], [], OfflineTraining.System, ROM.V);

end

%% ========================================================================
% Run POD-DEIM on nonlinear Snapshots
Fint_Reader  = HDF5_DenseMultiCVector(OfflineTraining.System.h5_filename, OfflineTraining.System.InternalForces.h5_section);
S_int        = Fint_Reader.readValues();
[Phi_int,    ~,    Sigma_int]  = VPOD_basis_computation(S_int, [], tol_POD_Fint, 1);
[IDEIM_int, P_int] = DEIM( Phi_int );

Fext_Reader  = HDF5_DenseMultiCVector(OfflineTraining.System.h5_filename, OfflineTraining.System.ExternalForces.h5_section);
S_ext        = Fext_Reader.readValues();
[Phi_ext,    ~,    Sigma_ext]  = VPOD_basis_computation(S_ext, [], tol_POD_Fext, 1);
[IDEIM_ext, P_ext] = DEIM( Phi_ext );

figure
subplot(1,2,1)
loglog(Sigma_int./Sigma_int(1),'-or');
hold on
cut_line_x = [1 : size(Phi_int,2)];
cut_line_y = Sigma_int(size(Phi_int,2))./Sigma_int(1) * ones(1,size(Phi_int,2));
loglog(cut_line_x, cut_line_y,'--k');
grid on
xlim([1  size(S_int,2)])
title('Fint snapshots spectrum')

subplot(1,2,2)
loglog(Sigma_ext./Sigma_ext(1),'-or');
hold on
cut_line_x = [1 : size(Phi_ext,2)];
cut_line_y = Phi_ext(size(Phi_ext,2))./Sigma_ext(1) * ones(1,size(Phi_ext,2));
loglog(cut_line_x, cut_line_y,'--k');
grid on
xlim([1  size(S_ext,2)])
title('Fext snapshots spectrum')

clear S_ext S_int;

ROM.Phi_int_IDEIM = Phi_int(IDEIM_int, :);
ROM.Phi_ext_IDEIM = Phi_ext(IDEIM_ext, :);
ROM.IDEIM_ext = IDEIM_ext;
ROM.IDEIM_in  = IDEIM_int;

ROM.LeftProjection_int = ( ROM.V' * Phi_int ) / ( ROM.Phi_int_IDEIM );
ROM.LeftProjection_ext = ( ROM.V' * Phi_ext ) / ( ROM.Phi_ext_IDEIM );

DATA       = CSM_read_DataFile('datafileR', dim, mu_bar);

% Set quad_order
if dim == 2
    quad_order       = 4;
elseif dim == 3
    quad_order       = 5;
end

[ MESH ] = buildMESH( dim, elements, vertices, boundaries, fem, quad_order, DATA, 'CSM' );

RedMeshObject =  ReducedMesh( MESH, fem, 'CSM' );
RedMeshObject.AppendInternalDoFs( IDEIM_int );
RedMeshObject.AppendInternalDoFs( IDEIM_ext );
RedMeshObject.Build( DATA );
RedMeshObject.ExportToVtk( 'FiguresTest/', 'ShearCube');

ROM.Red_Mesh = RedMeshObject.M_Red_Mesh;
            
%% ========================================================================
% Solve POD-DEIM ROM

%testing Parameters
mu_test_Dimension   = 5; % number of samples
mu_cube             = lhsdesign(mu_test_Dimension,P); % generate normalized design
Testing_Parameters  = bsxfun(@plus,mu_min,bsxfun(@times,mu_cube,(mu_max-mu_min)));


for i = 1 : mu_test_Dimension

    tmp_time = tic;
    U_ROM = CSM_PODDEIM_Solver(dim, elements, vertices, boundaries, fem, 'datafileR', ...
               Testing_Parameters(i,:), ['FiguresTest/TestingTierIII_', num2str(i)], [], ROM);
    time_ROM(i) = toc(tmp_time);
       
    tmp_time = tic;
    U_FEM = CSM_Solver(dim, elements, vertices, boundaries, fem, 'datafileR', ...
               Testing_Parameters(i,:));
    time_FOM(i) = toc(tmp_time);
    
    Error_Training(i) = norm(U_FEM(MESH.internal_dof) - U_ROM(MESH.internal_dof)) / norm( U_FEM(MESH.internal_dof) );      

end

fprintf('\n\n*** Average Relative Error on the displacement = %2.3f%% \n', mean( Error_Training ) * 100 );
fprintf('\n*** Average FOM Time-to-Solution = %1.2e s\n', mean(time_FOM) );
fprintf('\n*** Average ROM Time-to-Solution = %1.2e s\n\n', mean(time_ROM) );
%% ========================================================================

close all;

[~,~,~] = rmdir('SnapshotsTest', 's');
[~,~,~] = rmdir('FiguresTest', 's');

end
