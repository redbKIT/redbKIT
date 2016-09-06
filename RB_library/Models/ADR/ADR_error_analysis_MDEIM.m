function [ Error ] = ADR_error_analysis_MDEIM(FOM, ROM, N_samples, datafile, elements, vertices, boundaries, ROM_SEMMT)

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 


%% generate mu_sample
reset_randomSeed;
fprintf('Ntest = %d \n',N_samples);
for i = 1 : FOM.P
    sample(:,i) = rand(N_samples,1)*(FOM.mu_max(i)-FOM.mu_min(i))+FOM.mu_min(i);
end

N_vec = [1:2:ROM.N ROM.N];

delta_N       =   zeros(length(N_vec), N_samples);
delta_N_rel   =   zeros(length(N_vec), N_samples);
X_error       =   zeros(length(N_vec), N_samples);
X_error_rel   =   zeros(length(N_vec), N_samples);

Xnorm = FOM.Xnorm;


for k = 1 : N_samples
    
      mu = sample(k,:);
      fprintf('k_sample = %d \n',k);     
      
      def_nodes      =  MoveMesh(FOM.MESH.nodes, mu(:,FOM.DATA.shape_param), ROM_SEMMT);

      [uh] = Elliptic_Solver(FOM.MESH.dim,  elements, def_nodes(:,1:FOM.MESH.numVertices), boundaries, FOM.FE_SPACE.fem, datafile, mu);
      uh   = uh(FOM.MESH.internal_dof);
      
      %% RB solution and a posteriori error bounds      
      for n = 1:length(N_vec)
                        
            [~, uNh]        =  solve_RBsystem(ROM, mu, N_vec(n), [], [], def_nodes, ROM_SEMMT);
      
            delta_N(n,k)    =  1;%error_estimate(ROM, uN, mu);
            
            diff            =  uh - uNh(ROM.MESH.internal_dof);
            X_error(n,k)    =  sqrt(abs(diff'* (Xnorm * diff))); % absolute error
      end
      
      normU                  =  sqrt(abs(uh'* (FOM.Xnorm * uh)));
      X_error_rel(:,k)       =  X_error(:,k)./normU; % relative error
      delta_N_rel(:,k)       =  delta_N(:,k)./normU;
end

%% average and max error over mu_sample
X_error_av        =  mean(X_error,2);
X_error_rel_av    =  mean(X_error_rel,2);
delta_U_av        =  mean(delta_N,2);
delta_U_rel_av    =  mean(delta_N_rel,2);

X_error_max       =  max(X_error,[],2);
X_error_rel_max   =  max(X_error_rel,[],2);
delta_U_max       =  max(delta_N,[],2);
delta_U_rel_max   =  max(delta_N_rel,[],2);


[~,~,~] = mkdir('Figures');

%% plot absolute true errors and bounds
handle = figure;
semilogy(N_vec,X_error_av,'-ob','linewidth',2);
hold on
semilogy(N_vec,X_error_max,'-sb','linewidth',1);
grid on
title(['Absolute error vs N.  (n_{test}=',num2str(N_samples),')']);
legend('Average error','Max error');
xlabel('N');
set(findall(handle,'-property','FontSize'),'FontSize',14)

saveas(handle,'Figures/absolute_error_analysis','fig');

%% plot relative true errors and bounds
handle = figure;
semilogy(N_vec,X_error_rel_av,'-ob','linewidth',2);
hold on
semilogy(N_vec,X_error_rel_max,'-sb','linewidth',1);
grid on
title(['Relative error vs N.  (n_{test}=',num2str(N_samples),')']);
legend('Average error','Max error');
xlabel('N');
set(findall(handle,'-property','FontSize'),'FontSize',14)

saveas(handle,'Figures/relative_error_analysis','fig');


Error.sample                =  sample;

Error.X_error_av        =  X_error_av;
Error.X_error_rel_av    =  X_error_rel_av;
Error.X_error_max       =  X_error_max;
Error.Delta_N           =  delta_N;


end
