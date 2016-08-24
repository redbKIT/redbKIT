%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>   

clc
clear all
close all

[~,~,~] = mkdir('Snapshots');
[~,~,~] = mkdir('Figures');

delete('Snapshots/SystemSnapshots.h5')
delete('Snapshots/DisplacementSnapshots.h5')

%% ========================================================================
% DATA
dim      =  2;
fem      =  'P2';

mu_bar = [];

OfflineTraining.Solution.h5_filename       = 'Snapshots/DisplacementSnapshots.h5';
OfflineTraining.Solution.h5_section        = 'Displacement';
OfflineTraining.Solution.SamplingFrequency = 1;

OfflineTraining.System.h5_filename                      = 'Snapshots/SystemSnapshots.h5';
OfflineTraining.System.InternalForces.h5_section        = 'F_int';
OfflineTraining.System.InternalForces.SamplingFrequency = 1;
OfflineTraining.System.ExternalForces.h5_section        = 'F_ext';
OfflineTraining.System.ExternalForces.SamplingFrequency = 1;

tol_POD_U    = 1e-3;
tol_POD_Fext = 1e-3;
tol_POD_Fint = 1e-3;
%% ========================================================================
% Solve High-Fidelity Models and Collect solution Snapshots

% load P1 mesh
[vertices, boundaries, elements] = msh_to_Mmesh('mesh/Turek_mesh_Solid', dim);

[U, FE_SPACE, MESH, DATA] = CSMt_Solver(dim, elements, vertices, boundaries, fem, 'CSM3_data', ...
    [], [], true, OfflineTraining.Solution);

% PostProcessing
indexA = find(ismember(MESH.vertices(1:2,:)',[0.6 0.2],'rows'));

load ReferenceValues/csm3_l4_t0p005.point;

t = DATA.time.t0:DATA.time.dt:DATA.time.tf;

handle = figure;
subplot(2,1,1)
plot(t,U(indexA,1:length(t)))
hold on
plot(csm3_l4_t0p005(:,1),csm3_l4_t0p005(:,11),'--r')
legend('Test', 'Reference value')
xlabel('time [s]')
ylabel('x-Displacement [m]')
grid on


subplot(2,1,2)
plot(t,U(indexA+FE_SPACE.numDofScalar,1:length(t)))
hold on
plot(csm3_l4_t0p005(:,1),csm3_l4_t0p005(:,12),'--r')
legend('Test', 'Reference value')
xlabel('time [s]')
ylabel('y-Displacement [m]')
hold on
grid on

saveas(handle,'Figures/FOM_Adisplacement','epsc');
saveas(handle,'Figures/FOM_Adisplacement','fig');

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

[U_POD, Mass] = CSMt_POD_Solver(dim, elements, vertices, boundaries, fem, 'CSM3_data', ...
    [], [], true, OfflineTraining.System, ROM.V);

% PostProcessing
indexA = find(ismember(MESH.vertices(1:2,:)',[0.6 0.2],'rows'));

load ReferenceValues/csm3_l4_t0p005.point;

t = DATA.time.t0:DATA.time.dt:DATA.time.tf;

handle = figure;
subplot(2,1,1)
plot(t,U_POD(indexA,1:length(t)))
hold on
plot(csm3_l4_t0p005(:,1),csm3_l4_t0p005(:,11),'--r')
legend('Test', 'Reference value')
xlabel('time [s]')
ylabel('x-Displacement [m]')
grid on


subplot(2,1,2)
plot(t,U_POD(indexA+FE_SPACE.numDofScalar,1:length(t)))
hold on
plot(csm3_l4_t0p005(:,1),csm3_l4_t0p005(:,12),'--r')
legend('Test', 'Reference value')
xlabel('time [s]')
ylabel('y-Displacement [m]')
hold on
grid on

saveas(handle,'Figures/PODG_ROM_Adisplacement','epsc');
saveas(handle,'Figures/PODG_ROM_Adisplacement','fig');

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

ROM.M = ROM.V' * (Mass(MESH.internal_dof, MESH.internal_dof) * ROM.V);

DATA       = CSM_read_DataFile('CSM3_data', dim, mu_bar);

RedMeshObject =  ReducedMesh( MESH, fem, 'CSM' );
RedMeshObject.AppendInternalDoFs( IDEIM_int );
RedMeshObject.AppendInternalDoFs( IDEIM_ext );
RedMeshObject.Build( DATA );
RedMeshObject.ExportToVtk( 'Figures/', 'Turek2D');

ROM.Red_Mesh = RedMeshObject.M_Red_Mesh;


%% ========================================================================
% Solve POD-DEIM ROM

[U_PODDEIM, ~, ~, DATA_R] = CSMt_PODDEIM_Solver(dim, elements, vertices, boundaries, fem, 'CSM3_data', ...
    [], [], true, ROM);

% PostProcessing
indexA = find(ismember(MESH.vertices(1:2,:)',[0.6 0.2],'rows'));

load ReferenceValues/csm3_l4_t0p005.point;

t = DATA_R.time.t0:DATA_R.time.dt:DATA_R.time.tf;

handle = figure;
subplot(2,1,1)
plot(t,U_PODDEIM(indexA,1:length(t)))
hold on
plot(csm3_l4_t0p005(:,1),csm3_l4_t0p005(:,11),'--r')
legend('Test', 'Reference value')
xlabel('time [s]')
ylabel('x-Displacement [m]')
grid on


subplot(2,1,2)
plot(t,U_PODDEIM(indexA+FE_SPACE.numDofScalar,1:length(t)))
hold on
plot(csm3_l4_t0p005(:,1),csm3_l4_t0p005(:,12),'--r')
legend('Test', 'Reference value')
xlabel('time [s]')
ylabel('y-Displacement [m]')
hold on
grid on

saveas(handle,'Figures/PODDEIM_ROM_Adisplacement','epsc');
saveas(handle,'Figures/PODDEIM_ROM_Adisplacement','fig');

%% ========================================================================
