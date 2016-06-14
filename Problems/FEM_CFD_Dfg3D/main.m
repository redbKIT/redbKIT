%   Solution of a 3D Laminar Flow around a Cylinder benchmark problem.
%
%   Settings of the problem and reference values for the solution can be
%   found at
%     http://www.featflow.de/en/benchmarks/cfdbenchmarking/flow/dfg_flow3d.html
%
%   See also the following paper:
%   Solutions of 3D Navier-Stokes benchmark problems with adaptive finite elements
%   by M. Braack, T. Richter. Available at 
%     http://numerik.iwr.uni-heidelberg.de/Paper/braackrichter-200410.pdf
%
%   Both Reynolds = 20 and Reynolds = 100 cases are provided.

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

clc
clear all

test_case = 'Re100';

[~,~,~]  = mkdir('Results');
dim      =  3;

switch test_case
    
    case 'Re20'
        
        % load P1 mesh - several refinement levels are available
        [vertices, boundaries, elements] = msh_to_Mmesh('mesh/fluid_L1', dim);
        
        % specify FEM approximation
        fem        = {'P2', 'P1'};
        
        [U, MESH, DATA] = NSt_Solver(dim, elements, vertices, boundaries, fem, 'datafileRe20');
        
        % Postprocessing
        RefDrag = 6.185331;
        RefLift = 9.40135e-3;
        
        [t, Drag, Lift] = importAerodynamicForces( DATA.Output.DragLift.filename );
        
        fprintf('\nredbKIT RESULT:\n Drag = %1.3f, Lift = %1.5e\n',  Drag(end), Lift(end))
        fprintf('\nREFERENCE VALUES:\n Drag = %1.3f, Lift = %1.5e\n',   RefDrag,  RefLift)
        
        handle = figure;
        subplot(2,1,1)
        plot(t,Drag,'-b','LineWidth',2)
        hold on
        plot(t,Drag*0+RefDrag,'--k','LineWidth',2)
        legend('Drag Coefficient','Reference Value')
        xlabel('Time [s]')
        grid on
        
        subplot(2,1,2)
        plot(t,Lift,'-r','LineWidth',2)
        hold on
        plot(t,Lift*0+RefLift,'--k','LineWidth',2)
        legend('Lift Coefficient','Reference Value')
        xlabel('Time [s]')
        grid on
        
        saveas(handle,'Results/Re20_DragLift','epsc');
        saveas(handle,'Results/Re20_DragLift','fig');
        
    case 'Re100'
        
        % load P1 mesh - several refinement levels are available
        [vertices, boundaries, elements] = msh_to_Mmesh('mesh/fluid_L1', dim);
        
        % specify FEM approximation
        fem        = {'P1', 'P1'};
        
        [U, MESH, DATA] = NSt_Solver(dim, elements, vertices, boundaries, fem, 'datafileRe100');
        
        % Postprocessing
        [t, Drag, Lift, ZForce]     = importAerodynamicForces( DATA.Output.DragLift.filename );
        
        [t_ref,DragR,LiftR,ZForceR] = importBenchValues( 'ReferenceValues/BenchValues.txt');
                        
        handle = figure;
        subplot(3,1,1)
        plot(t,Drag,'-b','LineWidth',2)
        hold on
        plot(t_ref,DragR,'--k','LineWidth',2)
        legend('Drag Coefficient','Reference Value')
        xlabel('Time [s]')
        %xlim([3 max(t_ref)])
        grid on
        
        subplot(3,1,2)
        plot(t,Lift,'-r','LineWidth',2)
        hold on
        plot(t_ref,LiftR,'--k','LineWidth',2)
        legend('Lift Coefficient','Reference Value')
        xlabel('Time [s]')
        %xlim([5 max(t_ref)])
        grid on
        
        subplot(3,1,3)
        plot(t,ZForce,'-g','LineWidth',2)
        hold on
        plot(t_ref,ZForceR,'--k','LineWidth',2)
        legend('Z-Force','Reference Value')
        xlabel('Time [s]')
        %xlim([5 max(t_ref)])
        grid on
        
        saveas(handle,'Results/Re100_DragLift','epsc');
        saveas(handle,'Results/Re100_DragLift','fig');
        
end