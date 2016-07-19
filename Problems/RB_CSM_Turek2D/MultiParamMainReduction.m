%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>   

clc
clear all
close all

[~,~,~] = mkdir('Snapshots');
[~,~,~] = mkdir('Figures');

delete('Snapshots/SystemSnapshotsMP.h5')
delete('Snapshots/DisplacementSnapshotsMP.h5')

%% ========================================================================
% DATA
dim      =  2;
fem      =  'P1';

% Parameters range
P = 3;% Yuoung modulus, Poisson Coefficient, External load
mu_min = [ 0.5];
mu_max = [ 4];
mu_bar = [ 1]; 
    
% Training Parameters
mu_train_Dimension   = 3; % number of samples
Training_Parameters  = [0.5 2 4];

h5_filename_Sol       = 'Snapshots/DisplacementSnapshotsMP';
h5_filename_Sys       = 'Snapshots/SystemSnapshotsMP';

OfflineTraining.Solution.h5_section        = 'Displacement';
OfflineTraining.Solution.SamplingFrequency = 1;

OfflineTraining.System.InternalForces.h5_section        = 'F_int';
OfflineTraining.System.InternalForces.SamplingFrequency = 1;
OfflineTraining.System.ExternalForces.h5_section        = 'F_ext';
OfflineTraining.System.ExternalForces.SamplingFrequency = 1;

tol_POD_U    = 1e-3;
tol_POD_Fext = 1e-3;
tol_POD_Fint = 1e-3;

tol_POD_U_local    = 1e-3;
tol_POD_Fext_local = 1e-3;
tol_POD_Fint_local = 1e-3;
%% ========================================================================
% Solve High-Fidelity Models and Collect solution Snapshots

% load P1 mesh
[vertices, boundaries, elements] = msh_to_Mmesh('mesh/Turek_mesh_Solid', dim);

for i = 1 : mu_train_Dimension
    
    OfflineTraining.Solution.h5_filename       = [h5_filename_Sol, num2str(i), '.h5'];

    [U, FE_SPACE, MESH, DATA] = CSMt_Solver(dim, elements, vertices, boundaries, fem, 'datafileMP', ...
        Training_Parameters(i), [], true, OfflineTraining.Solution);
    
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
    
    saveas(handle,['Figures/FOM_Adisplacement',num2str(i)],'epsc');
    saveas(handle,['Figures/FOM_Adisplacement',num2str(i)],'fig');
    
end

%% ========================================================================
% Generate POD basis for the solution

V = [];

for i = 1 : mu_train_Dimension
    DispSnap_HDF5  = HDF5_DenseMultiCVector([h5_filename_Sol, num2str(i), '.h5'], OfflineTraining.Solution.h5_section);
    S_u            = DispSnap_HDF5.readValues();
    V_local        = VPOD_basis_computation(S_u, [], tol_POD_U_local, 1);
    V              = VPOD_basis_computation([V V_local], [], tol_POD_U, 1);
end

ROM.V = V;

%% ========================================================================
% Solve POD-Galerkin ROM and Collect System Snapshots

for i = 1 : mu_train_Dimension
    
    OfflineTraining.System.h5_filename       = [h5_filename_Sys, num2str(i), '.h5'];
    
    [U_POD, Mass] = CSMt_POD_Solver(dim, elements, vertices, boundaries, fem, 'datafileMP', ...
        Training_Parameters(i), [], true, OfflineTraining.System, ROM.V);
    
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
    
    saveas(handle,['Figures/PODG_ROM_Adisplacement',num2str(i)],'epsc');
    saveas(handle,['Figures/PODG_ROM_Adisplacement',num2str(i)],'fig');
    
end

%% ========================================================================
% Run POD-DEIM on nonlinear Snapshots

Phi_int = [];
for i = 1 : mu_train_Dimension
    
    Fint_Reader  = HDF5_DenseMultiCVector([h5_filename_Sys, num2str(i), '.h5'], OfflineTraining.System.InternalForces.h5_section);
    S_int        = Fint_Reader.readValues();
    Phi_int_loc  = VPOD_basis_computation(S_int, [], tol_POD_Fint_local, 1);
    Phi_int      = VPOD_basis_computation([Phi_int Phi_int_loc], [], tol_POD_Fint, 1);

end
clear S_int;
[IDEIM_int, P_int] = DEIM( Phi_int );

Phi_ext = [];
for i = 1 : mu_train_Dimension
    
    Fext_Reader  = HDF5_DenseMultiCVector([h5_filename_Sys, num2str(i), '.h5'], OfflineTraining.System.ExternalForces.h5_section);
    S_ext        = Fext_Reader.readValues();
    Phi_ext_loc  = VPOD_basis_computation(S_ext, [], tol_POD_Fext_local, 1);
    Phi_ext      = VPOD_basis_computation([Phi_ext Phi_ext_loc], [], tol_POD_Fext, 1);
    
end
clear S_ext;

[IDEIM_ext, P_ext] = DEIM( Phi_ext );

ROM.Phi_int_IDEIM = Phi_int(IDEIM_int, :);
ROM.Phi_ext_IDEIM = Phi_ext(IDEIM_ext, :);
ROM.IDEIM_ext = IDEIM_ext;
ROM.IDEIM_in  = IDEIM_int;

ROM.LeftProjection_int = ( ROM.V' * Phi_int ) / ( ROM.Phi_int_IDEIM );
ROM.LeftProjection_ext = ( ROM.V' * Phi_ext ) / ( ROM.Phi_ext_IDEIM );

ROM.M = ROM.V' * (Mass(MESH.internal_dof, MESH.internal_dof) * ROM.V);

DATA       = CSM_read_DataFile('datafileMP', dim, mu_bar);

RedMeshObject =  ReducedMesh( MESH, fem, 'CSM' );
RedMeshObject.AppendInternalDoFs( IDEIM_int );
RedMeshObject.AppendInternalDoFs( IDEIM_ext );
RedMeshObject.Build( DATA );
RedMeshObject.ExportToVtk( 'FiguresMP/', 'Test');

ROM.Red_Mesh = RedMeshObject.M_Red_Mesh;


%% ========================================================================
% Solve POD-DEIM ROM

[U_PODDEIM, ~, ~, DATA_R] = CSMt_PODDEIM_Solver(dim, elements, vertices, boundaries, fem, 'datafileMP', ...
    2, [], true, ROM);

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

saveas(handle,'Figures/PODDEIM_ROM_Adisplacement1','epsc');
saveas(handle,'Figures/PODDEIM_ROM_Adisplacement1','fig');

%% ========================================================================
