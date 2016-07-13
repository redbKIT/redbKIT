%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri@epfl.ch>

clc
clear all

dim      =  2;
fem      =  'P2';

[~,~,~] = mkdir('Figures');

%% Load P1 mesh
[vertices, boundaries, elements] = msh_to_Mmesh('../mesh/Turek_mesh_Solid', dim);

%% CSM1
fprintf('\n------------ TEST CSM 1 --------------\n')

[U, FE_SPACE, MESH, DATA] = CSM_Solver(dim, elements, vertices, boundaries, fem, 'CSM1_data', [], 'Figures/Sol_CSM1');

% PostProcessing
indexA = find(ismember(MESH.vertices(1:2,:)',[0.6 0.2],'rows'));

fprintf('\nCSM1 RESULT:\nA-Point: x-Displacement = %1.3e m, y-Displacement = %1.3e m\n', U(indexA), U(indexA+FE_SPACE.numDofScalar))
fprintf('\nREFERENCE VALUES:\nA-Point: x-Displacement = %1.3e m, y-Displacement = %1.3e m\n', -7.187*1e-3, -66.10*1e-3)


%% CSM2
fprintf('\n------------ TEST CSM 2 --------------\n')

[U, FE_SPACE, MESH, DATA] = CSM_Solver(dim, elements, vertices, boundaries, fem, 'CSM2_data', [], 'Figures/Sol_CSM2');

% PostProcessing
indexA = find(ismember(MESH.vertices(1:2,:)',[0.6 0.2],'rows'));

fprintf('\nCSM2 RESULT:\nA-Point: x-Displacement = %1.3e m, y-Displacement = %1.3e m\n', U(indexA), U(indexA+FE_SPACE.numDofScalar))
fprintf('\nREFERENCE VALUES:\nA-Point: x-Displacement = %1.3e m, y-Displacement = %1.3e m\n', -0.4690*1e-3, -16.97*1e-3)


%% CSM3
fprintf('\n------------ TEST CSM 3 --------------\n')

[U, FE_SPACE, MESH, DATA] = CSMt_Solver(dim, elements, vertices, boundaries, fem, 'CSM3_data', [], [], true);

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


saveas(handle,'Figures/CSM3_Adisplacement','epsc');
saveas(handle,'Figures/CSM3_Adisplacement','fig');

