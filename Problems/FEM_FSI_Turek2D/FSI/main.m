%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri@epfl.ch>

clear all
clc

dim      =  2;

[~,~,~] = mkdir('Figures');
[~,~,~] = mkdir('Results');

%% FSI1
fprintf('\n------------ TEST FSI 1 --------------\n')

[mshS.vertices, mshS.boundaries, mshS.elements, mshS.rings] = msh_to_Mmesh('../mesh/Turek_mesh_Solid', dim);
[mshF.vertices, mshF.boundaries, mshF.elements, mshF.rings] = msh_to_Mmesh('../mesh/Turek_mesh_Fluid', dim);

[U1, MESH, DATA] = FSIt_Solver(dim, mshF, mshS, {'P2','P1'}, 'P2', 'FSI_TurekC1_NS_data', 'FSI_TurekC1_CSM_data', [], 'Figures/Turek1_');

% PostProcessing
[t, Drag, Lift] = importAerodynamicForces( DATA.Fluid.Output.DragLift.filename );

fprintf('\nFSI-1 RESULT:\n Drag = %2.3f, Lift = %1.3f\n',        Drag(end),    Lift(end))
fprintf('\nREFERENCE VALUES:\n Drag = %2.3f, Lift = %1.3f\n',    14.295,         0.7638)           
           
indexA = find(ismember(MESH.Solid.nodes(1:2,:)',[0.6 0.2],'rows'));

fprintf('\nFSI-1 RESULT:\nA-Point: x-Displacement = %1.3e m, y-Displacement = %1.3e m\n', ...
    U1(dim*MESH.Fluid.numNodes+MESH.Fluid.numVertices+indexA,end), U1(dim*MESH.Fluid.numNodes+MESH.Fluid.numVertices+indexA+MESH.Solid.numNodes,end))
fprintf('\nREFERENCE VALUES:\nA-Point: x-Displacement = %1.3e m, y-Displacement = %1.3e m\n', 0.0227*1e-3, 0.8209*1e-3)


%% FSI2
fprintf('\n------------ TEST FSI 2 --------------\n')

[mshS.vertices, mshS.boundaries, mshS.elements, mshS.rings] = msh_to_Mmesh('../mesh/Turek_mesh_Solid', dim);
[mshF.vertices, mshF.boundaries, mshF.elements, mshF.rings] = msh_to_Mmesh('../mesh/Turek_mesh_Fluid', dim);

[U2, MESH, DATA] = FSIt_Solver(dim, mshF, mshS, {'P1','P1'}, 'P1', 'FSI_TurekC2_NS_data', 'FSI_TurekC2_CSM_data', [], 'Figures/Turek2_');

% PostProcessing
[t, Drag, Lift] = importAerodynamicForces( DATA.Fluid.Output.DragLift.filename );

indexA = find(ismember(MESH.Solid.nodes(1:2,:)',[0.6 0.2],'rows'));

load ReferenceValues/ref_fsi2;
t_ref   = ref_fsi2(:,1);

handle = figure;
subplot(2,2,1)
plot(t,Drag,'-b','LineWidth',2)
hold on
plot(t_ref,ref_fsi2(:,5)+ref_fsi2(:,7),'--k','LineWidth',2)
legend('Drag Coefficient','Reference Value')
xlabel('Time [s]')
grid on

subplot(2,2,2)
plot(t,Lift,'-r','LineWidth',2)
hold on
plot(t_ref,ref_fsi2(:,6)+ref_fsi2(:,8),'--k','LineWidth',2)
legend('Lift Coefficient','Reference Value')
xlabel('Time [s]')
grid on


subplot(2,2,3)
plot(t,U2(dim*MESH.Fluid.numNodes+MESH.Fluid.numVertices+indexA,:),'-m','LineWidth',2)
hold on
plot(t_ref,ref_fsi2(:,11),'--r')
legend('Test', 'Reference value')
xlabel('Time [s]')
ylabel('x-Displacement [m]')
grid on


subplot(2,2,4)
plot(t,U2(dim*MESH.Fluid.numNodes+MESH.Fluid.numVertices+indexA+MESH.Solid.numNodes,:),'-m','LineWidth',2)
hold on
plot(t_ref,ref_fsi2(:,12),'--r')
legend('Test', 'Reference value')
xlabel('Time [s]')
ylabel('y-Displacement [m]')
grid on

saveas(handle,'Results/FSI2_DragLift','epsc');
saveas(handle,'Results/FSI2_DragLift','fig');
