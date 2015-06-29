function [ Error ] = error_analysis_mu(FOM, ROM, N_samples, N, mu_index, mu_fixed, value_mu_fixed)
%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

n_N = length(N);

%% generate mu_sample
reset_randomSeed;
ii_fix = 1;
for i = mu_fixed
      if nargin > 6
            val     =  value_mu_fixed(ii_fix);
            ii_fix  =  ii_fix + 1;
      else
            val     = rand(1,1)*(FOM.mu_max(i)-FOM.mu_min(i))+FOM.mu_min(i);
      end
      sample(:,i) = val*ones(N_samples,1);
end

for i = mu_index
      sample(:,i) = linspace(FOM.mu_min(i), FOM.mu_max(i), N_samples);
end

if isfield(ROM, 'Greedy_samples')
    sample = [sample; ROM.Greedy_samples];
end

N_samples     = size(sample,1);
fprintf('Ntest = %d \n',N_samples);

delta_N       =   zeros(n_N, N_samples);
delta_N_rel   =   zeros(n_N, N_samples);
X_error       =   zeros(n_N, N_samples);
X_error_rel   =   zeros(n_N, N_samples);


parfor k = 1 : N_samples
    
      mu = sample(k,:);
      fprintf('k_sample = %d \n',k);     
      
      uh = solve_HFsystem(FOM, mu);
      uh = uh(FOM.MESH.internal_dof);
            
      %% RB solution and a posteriori error bound     
      for n = 1 : n_N
          [uN, uNh]        =  solve_RBsystem(ROM, mu, N(n));
          delta_N(n, k)    =  error_estimate(ROM, uN, mu);
          
          diff             =  uh - uNh(FOM.MESH.internal_dof);
          X_error(n, k)    =  sqrt(abs(diff'* (FOM.Xnorm * diff))); % absolute error
          
          normU                =  sqrt(uh'* (FOM.Xnorm * uh));
          X_error_rel(n, k)    =  X_error(n, k)/normU; % relative error
          delta_N_rel(n, k)    =  delta_N(n, k)/normU;
      end
end

%% sorting for plot
[s, Is]      =  sort(sample(:,mu_index(1)));
sample       =  sample(Is,:);
X_error      =  X_error(:, Is);
delta_N      =  delta_N(:, Is);
X_error_rel  =  X_error_rel(:, Is);
delta_N_rel  =  delta_N_rel(:, Is);

[~,~,~] = mkdir('Figures');

%% plot relative true error and bound
for n = 1 : n_N
    h1=figure;
    for k = 1 : size(sample,2)
        subplot(ceil(size(sample,2)/2),2,k);
        semilogy(sample(:,k),X_error_rel(n, :),'-b', 'linewidth',2);
        hold on
        semilogy(sample(:,k), delta_N_rel(n, :), '-k', 'linewidth',2);
        grid on
        legend('Error','\Delta_N');
        set(findall(h1,'-property','FontSize'),'FontSize',14)
    end
    
    saveas(h1,['Figures/mu_error_analysis_N',num2str(N(n))],'fig');
    saveas(h1,['Figures/mu_error_analysis_N',num2str(N(n))],'epsc')
end

Error.sample        =  sample;
Error.X_error       =  X_error;
Error.Delta_N       =  delta_N;


end
