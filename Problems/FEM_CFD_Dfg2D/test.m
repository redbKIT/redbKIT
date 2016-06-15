function [] = test( velocity_FE_space, nonlinear_treatment, Stabilization )
%   Solution of a 2D Laminar Flow around a Cylinder benchmark problem.
%   See also main.m

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>


[~,~,~]  = mkdir('ResultsTest');
dim      =  2;

% load P1 mesh
[vertices, boundaries, elements] = msh_to_Mmesh('mesh/obstacleDFG_Coarse', dim);

% P2-P1 approximation
fem        = {velocity_FE_space, 'P1'};

ForcesFilename = ['ResultsTest/AerodynamicForces', velocity_FE_space, nonlinear_treatment, Stabilization,'.txt'];

[U, MESH, DATA] = NSt_Solver(dim, elements, vertices, boundaries, fem, 'datafileTest', {nonlinear_treatment, Stabilization, ForcesFilename});

% Postprocessing
[t, Drag, Lift] = importAerodynamicForces( DATA.Output.DragLift.filename );

load ReferenceValues_Re100Fixed/bdforces_lv6;

t_ref  = bdforces_lv6(:,2);

handle = figure;
subplot(2,1,1)
plot(t,Drag,'-b','LineWidth',2)
hold on
plot(t_ref,bdforces_lv6(:,4),'--k','LineWidth',2)
legend('Drag Coefficient','Reference Value')
xlabel('Time [s]')
xlim([0 4])
grid on

subplot(2,1,2)
plot(t,Lift,'-r','LineWidth',2)
hold on
plot(t_ref,bdforces_lv6(:,5),'--k','LineWidth',2)
legend('Lift Coefficient','Reference Value')
xlabel('Time [s]')
xlim([0 4])
grid on

saveas(handle,['ResultsTest/',velocity_FE_space, nonlinear_treatment, Stabilization, '_DragLift'],'epsc');
saveas(handle,['ResultsTest/',velocity_FE_space, nonlinear_treatment, Stabilization, '_DragLift'],'fig');


end