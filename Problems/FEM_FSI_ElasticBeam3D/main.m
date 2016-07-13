clear all
clc

[~,~,~] = mkdir('Figures');
[~,~,~] = mkdir('Results');

dim      =  3;

%% Load meshes
[mshS.vertices, mshS.boundaries, mshS.elements, mshS.rings] = msh_to_Mmesh('mesh/mesh_Solid', dim);
[mshF.vertices, mshF.boundaries, mshF.elements, mshF.rings] = msh_to_Mmesh('mesh/mesh_Fluid', dim);

%% Solve
[U2, MESH, DATA] = FSIt_Solver(dim, mshF, mshS, {'P1','P1'}, 'P1', 'NS_data', 'CSM_data', [], 'Figures/Benchmark_ImplP1NH_Final_');

%% PostProcessing
[t, Drag, Lift] = importAerodynamicForces( DATA.Fluid.Output.DragLift.filename );

% indexA = find(ismember(MESH.Solid.nodes(1:3,:)',[9.5 6],'rows'));
% 
% handle = figure;
% subplot(3,1,1)
% plot(t,Drag,'-b','LineWidth',2)
% legend('Drag Coefficient')
% xlabel('Time [s]')
% grid on
% 
% subplot(3,1,2)
% plot(t,Lift,'-r','LineWidth',2)
% legend('Lift Coefficient')
% xlabel('Time [s]')
% grid on
% 
% 
% subplot(3,1,3)
% plot(t,U2(dim*MESH.Fluid.numNodes+MESH.Fluid.numVertices+indexA+MESH.Solid.numNodes,:),'-m','LineWidth',2)
% xlabel('Time [s]')
% ylabel('Tip y-Displacement [m]')
% grid on
% 
% saveas(handle,'Results/FSI_DragLiftDisp','epsc');
% saveas(handle,'Results/FSI_DragLiftDisp','fig');