function [ Error ] = EIM_error_analysis(FOM, ROM, N_samples, data_file, N_vec, varargin)
%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

if nargin < 5 || isempty(N_vec)
    N_vec         =   1 : 1 : ROM.N;
end

%% generate mu_sample
reset_randomSeed;
fprintf('Ntest = %d \n',N_samples);
for i = 1 : FOM.P
    sample(:,i) = rand(N_samples,1)*(FOM.mu_max(i)-FOM.mu_min(i))+FOM.mu_min(i);
end

delta_N       =   zeros(length(N_vec), N_samples);
delta_N_rel   =   zeros(length(N_vec), N_samples);
X_error       =   zeros(length(N_vec), N_samples);
X_error_rel   =   zeros(length(N_vec), N_samples);


for k = 1 : N_samples
    
      mu = sample(k,:);
      fprintf('k_sample = %d \n',k);     
      
      [uh] = Elliptic_Solver(FOM.MESH.dim, FOM.MESH.elements, FOM.MESH.vertices, FOM.MESH.boundaries, FOM.FE_SPACE.fem, data_file, mu);
      %uh = solve_HFsystem(FOM, mu);
      uh  = uh(FOM.MESH.internal_dof);
            
      %% RB solution and a posteriori error bounds      
      for n = 1 : length(N_vec)
            [uN, uNh]       =  solve_RBsystem(ROM, mu, N_vec(n), varargin{:});
            delta_N(n,k)    =  error_estimate(ROM, uN, mu, varargin{:});
            
            diff            =  uh - uNh(FOM.MESH.internal_dof);
            X_error(n,k)    =  sqrt(abs(diff'* (FOM.Xnorm * diff))); % absolute error
      end
      
      normU                  =  sqrt(uh'* (FOM.Xnorm * uh));
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

if nargin > 5
    addM = ['_', num2str(varargin{1})] ; 
else
    addM = [];
end

%% plot absolute true errors and bounds
n = N_vec;
handle = figure;
semilogy(n,X_error_av,'-ob',n,delta_U_av,'-ok','linewidth',2);
hold on
semilogy(n,X_error_max,'-sb',n,delta_U_max,'-sk','linewidth',1);
grid on
title(['Absolute a posteriori error bounds and computed errors vs N.  (n_{test}=',num2str(N_samples),')']);
legend('Average error','\Delta_N average','Max error','\Delta_N max');
xlabel('N');
set(findall(handle,'-property','FontSize'),'FontSize',14)

saveas(handle,['Figures/absolute_error_analysis',addM],'fig');

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

saveas(handle,['Figures/relative_error_analysis',addM],'fig');


Error.sample                =  sample;

Error.X_error_av        =  X_error_av;
Error.X_error_rel_av    =  X_error_rel_av;
Error.X_error_max       =  X_error_max;
Error.Delta_N           =  delta_N;
Error.X_error           =  X_error;

end
