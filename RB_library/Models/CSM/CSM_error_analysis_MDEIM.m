function [ Error ] = CSM_error_analysis_MDEIM(FOM, ROM, N_samples, datafile, elements, vertices, boundaries)
%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 


%% generate mu_sample
reset_randomSeed;
fprintf('Ntest = %d \n',N_samples);
for i = 1 : FOM.P
    sample(:,i) = rand(N_samples,1)*(FOM.mu_max(i)-FOM.mu_min(i))+FOM.mu_min(i);
end

delta_N       =   zeros(ROM.N, N_samples);
delta_N_rel   =   zeros(ROM.N, N_samples);
X_error       =   zeros(ROM.N, N_samples);
X_error_rel   =   zeros(ROM.N, N_samples);


for k = 1 : N_samples
    
      mu = sample(k,:);
      fprintf('k_sample = %d \n',k);     
      
      uh = CSM_Solver(FOM.MESH.dim,  elements, vertices, boundaries, FOM.FE_SPACE.fem, datafile, mu);
      uh = uh(FOM.MESH.internal_dof);
   
      %% RB solution and a posteriori error bounds      
      for n = 1 : ROM.N
            [uN, uNh]       =  solve_RBsystem(ROM, mu, n);
            delta_N(n,k)    =  error_estimate(ROM, uN, mu);%error_estimate_full(FOM, ROM.V(:,1:n)*uN, mu);%
            
            diff            =  uh - uNh(FOM.MESH.internal_dof);
            X_error(n,k)    =  sqrt(abs(diff'* (FOM.Xnorm * diff))); % absolute error
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
n = 1:ROM.N;
handle = figure;
semilogy(n,X_error_av,'-ob',n,delta_U_av,'-ok','linewidth',2);
hold on
semilogy(n,X_error_max,'-sb',n,delta_U_max,'-sk','linewidth',1);
grid on
title(['Absolute a posteriori error bounds and computed errors vs N.  (n_{test}=',num2str(N_samples),')']);
legend('Average error','\Delta_N average','Max error','\Delta_N max');
xlabel('N');
set(findall(handle,'-property','FontSize'),'FontSize',14)

saveas(handle,'Figures/absolute_error_analysis','fig');

%% plot relative true errors and bounds
handle = figure;
semilogy(n,X_error_rel_av,'-ob',n,delta_U_rel_av,'-ok','linewidth',2);
hold on
semilogy(n,X_error_rel_max,'-sb',n,delta_U_rel_max,'-sk','linewidth',1);
grid on
title(['Relative a posteriori error bounds and computed errors vs N.  (n_{test}=',num2str(N_samples),')']);
legend('Average error','\Delta_N average','Max error','\Delta_N max');
xlabel('N');
set(findall(handle,'-property','FontSize'),'FontSize',14)

saveas(handle,'Figures/relative_error_analysis','fig');


Error.sample                =  sample;

Error.X_error_av        =  X_error_av;
Error.X_error_rel_av    =  X_error_rel_av;
Error.X_error_max       =  X_error_max;
Error.Delta_N           =  delta_N;


end
