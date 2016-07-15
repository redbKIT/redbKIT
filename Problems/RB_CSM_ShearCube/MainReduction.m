%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>   

clc
clear all

[~,~,~] = mkdir('Snapshots');
[~,~,~] = mkdir('Figures');

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

OfflineTraining.Solution.h5_filename       = 'Snapshots/DisplacementSnapshots.h5';
OfflineTraining.Solution.h5_section        = 'Displacement';
%OfflineTraining.Solution.SamplingFrequency = 1;

OfflineTraining.System.h5_filename                      = 'Snapshots/SystemSnapshots.h5';
OfflineTraining.System.InternalForces.h5_section        = 'F_int';
%OfflineTraining.System.InternalForces.SamplingFrequency = 1;
OfflineTraining.System.ExternalForces.h5_section        = 'F_ext';
%OfflineTraining.System.ExternalForces.SamplingFrequency = 1;

tol_POD_U    = 1e-4;
tol_POD_Fext = 1e-4;
tol_POD_Fint = 1e-4;
%% ========================================================================
% Solve High-Fidelity Models and Collect solution Snapshots

% load P1 mesh
[vertices, boundaries, elements] = msh_to_Mmesh('Cube', dim);

for i = 1 : mu_train_Dimension

    CSM_Solver(dim, elements, vertices, boundaries, fem, 'datafileR', ...
               Training_Parameters(i,:), ['Figures/TrainingTierI_', num2str(i)], [], OfflineTraining.Solution);

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

%% ========================================================================
% Solve POD-Galerkin ROM and Collect System Snapshots

for i = 1 : mu_train_Dimension

    CSM_POD_Solver(dim, elements, vertices, boundaries, fem, 'datafileR', ...
               Training_Parameters(i,:), ['Figures/TrainingTierII_', num2str(i)], [], OfflineTraining.System, V);

end

%% ========================================================================
% Run POD-DEIM on nonlinear Snapshots
Fint_Reader  = HDF5_DenseMultiCVector(OfflineTraining.System.h5_filename, OfflineTraining.System.InternalForces.h5_section);
S_int        = Fint_Reader.readValues();
[Phi_int,    ~,    Sigma_int]  = VPOD_basis_computation(S_int, [], tol_POD_Fint, 1);

Fext_Reader  = HDF5_DenseMultiCVector(OfflineTraining.System.h5_filename, OfflineTraining.System.ExternalForces.h5_section);
S_ext        = Fext_Reader.readValues();
[Phi_ext,    ~,    Sigma_ext]  = VPOD_basis_computation(S_ext, [], tol_POD_Fext, 1);

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
cut_line_y = Phi_ext(size(V,2))./Sigma_ext(1) * ones(1,size(Phi_ext,2));
loglog(cut_line_x, cut_line_y,'--k');
grid on
xlim([1  size(S_ext,2)])
title('Fext snapshots spectrum')

%% ========================================================================
