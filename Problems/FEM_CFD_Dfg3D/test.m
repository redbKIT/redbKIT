function [] = test( velocity_FE_space, nonlinear_treatment, Stabilization )
%   Solution of a 3D Laminar Flow around a Cylinder benchmark problem.
%   See also main.m

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>


[~,~,~]  = mkdir('ResultsTest');
dim      = 3;

% load P1 mesh
[vertices, boundaries, elements] = msh_to_Mmesh('mesh/fluid_L0', dim);

% P2-P1 approximation
fem        = {velocity_FE_space, 'P1'};

ForcesFilename = ['ResultsTest/AerodynamicForces', velocity_FE_space, nonlinear_treatment, Stabilization,'.txt'];

[~, ~, DATA] = NSt_Solver(dim, elements, vertices, boundaries, fem, 'datafileTest', {nonlinear_treatment, Stabilization, ForcesFilename});

% Postprocessing
[t, Drag, Lift, ZForce] = importAerodynamicForces( DATA.Output.DragLift.filename );

[t_ref,DragR,LiftR,ZForceR] = importBenchValues( 'ReferenceValues/BenchValues.txt');

handle = figure;
subplot(3,1,1)
plot(t,Drag,'-ob','LineWidth',2)
hold on
plot(t_ref,DragR,'--k','LineWidth',2)
legend('Drag Coefficient','Reference Value')
xlabel('Time [s]')
xlim([0 1])
grid on

subplot(3,1,2)
plot(t,Lift,'-or','LineWidth',2)
hold on
plot(t_ref,LiftR,'--k','LineWidth',2)
legend('Lift Coefficient','Reference Value')
xlabel('Time [s]')
xlim([0 1])
grid on

subplot(3,1,3)
plot(t,ZForce,'-og','LineWidth',2)
hold on
plot(t_ref,ZForceR,'--k','LineWidth',2)
legend('Z-Force','Reference Value')
xlabel('Time [s]')
xlim([0 1])
grid on

saveas(handle,['ResultsTest/',velocity_FE_space, nonlinear_treatment, Stabilization, '_DragLift'],'epsc');
saveas(handle,['ResultsTest/',velocity_FE_space, nonlinear_treatment, Stabilization, '_DragLift'],'fig');


end