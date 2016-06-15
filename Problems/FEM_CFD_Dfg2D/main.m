%   Solution of a 2D Laminar Flow around a Cylinder benchmark problem.
% 
%   "Flow around a cylinder" is a well-known benchmark for the evaluation of 
%   numerical algorithms for incompressible Navier-Stokes equations in the laminar case. 
%   It was developed in 1995 as part of the research project 
%   "Flow simulation on high-performance computers" funded by the German 
%   Research Association (DFG).
%
%   Settings of the problem and reference values for the solution can be
%   found at
%     http://www.featflow.de/en/benchmarks/cfdbenchmarking/flow.html
%
%   See also the following paper:
%   Benchmark computations of laminar flow around cylinder,
%   by M. Shafer, S. Turek. Available at 
%     http://www.math.u-psud.fr/~maury/paps/NUM_NS.BenchmarkTurek.pdf
%
%   Both Reynolds = 20 and Reynolds = 100 (periodic and fixed time interval)
%   cases are provided.

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>


clc
clear all

test_case = 'Re100Periodic';

[~,~,~]  = mkdir('Results');
dim      =  2;

switch test_case
    
    case 'Re20'
        
        % load P1 mesh
        [vertices, boundaries, elements] = msh_to_Mmesh('mesh/obstacleDFG_Coarse', dim);
        
        % P2-P1 approximation
        fem        = {'P2', 'P1'};
        
        [U, MESH, DATA] = NSt_Solver(dim, elements, vertices, boundaries, fem, 'datafileRe20');
        
        % Postprocessing
        RefDrag = 5.579;
        RefLift = 0.0106;
        RefDP   = 0.11752;
        
        indexA1 = find(ismember(vertices(1:2,:)',[0.15 0.2],'rows'));
        indexA2 = find(ismember(vertices(1:2,:)',[0.25 0.2],'rows'));
        
        numNodes = MESH.numNodes;
        DeltaP   = U(dim*numNodes+indexA1) - U(dim*numNodes+indexA2);
        [t, Drag, Lift] = importAerodynamicForces( DATA.Output.DragLift.filename );
        
        fprintf('\nCFD1 RESULT:\n Drag = %1.3f, Lift = %1.4f, Delta P = %1.4f\n',  Drag(end), Lift(end), DeltaP(end))
        fprintf('\nREFERENCE VALUES:\n Drag = %1.3f, Lift = %1.4f, Delta P = %1.4f\n',   RefDrag,  RefLift, RefDP)
        
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
        
    case 'Re100Periodic'

        % load P1 mesh
        [vertices, boundaries, elements] = msh_to_Mmesh('mesh/obstacleDFG', dim);
        
        % P2-P1 approximation
        fem        = {'P2', 'P1'};
        
        [U, MESH, DATA] = NSt_Solver(dim, elements, vertices, boundaries, fem, 'datafileRe100Periodic');
        
        % Postprocessing
        [t, Drag, Lift] = importAerodynamicForces( DATA.Output.DragLift.filename );
        
        load ReferenceValues_Re100Periodic/bdforces_q2_lv6_dt4;
        
        t_ref  = bdforces_q2_lv6_dt4(:,2);
        
        handle = figure;
        subplot(2,1,1)
        plot(t,Drag,'-b','LineWidth',2)
        hold on
        plot(t_ref,bdforces_q2_lv6_dt4(:,4),'--k','LineWidth',2)
        legend('Drag Coefficient','Reference Value')
        xlabel('Time [s]')
        %xlim([3 max(t_ref)])
        grid on
        
        subplot(2,1,2)
        plot(t,Lift,'-r','LineWidth',2)
        hold on
        plot(t_ref,bdforces_q2_lv6_dt4(:,5),'--k','LineWidth',2)
        legend('Lift Coefficient','Reference Value')
        xlabel('Time [s]')
        %xlim([5 max(t_ref)])
        grid on
        
        saveas(handle,'Results/Re100Periodic_DragLift','epsc');
        saveas(handle,'Results/Re100Periodic_DragLift','fig');
        
        
    case 'Re100Fixed'
        
        % load P1 mesh
        [vertices, boundaries, elements] = msh_to_Mmesh('mesh/obstacleDFG', dim);
        
        % P2-P1 approximation
        fem        = {'P2', 'P1'};
        
        [U, MESH, DATA] = NSt_Solver(dim, elements, vertices, boundaries, fem, 'datafileRe100Fixed');
        
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
        grid on
        
        subplot(2,1,2)
        plot(t,Lift,'-r','LineWidth',2)
        hold on
        plot(t_ref,bdforces_lv6(:,5),'--k','LineWidth',2)
        legend('Lift Coefficient','Reference Value')
        xlabel('Time [s]')
        grid on
        
        saveas(handle,'Results/Re100Fixed_DragLift','epsc');
        saveas(handle,'Results/Re100Fixed_DragLift','fig');
        
end
