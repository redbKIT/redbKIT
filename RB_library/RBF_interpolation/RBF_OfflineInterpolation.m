function [FOM] = RBF_OfflineInterpolation(FOM)
%RBF_OFFLINEINTERPOLATION offline construction of the RBF interpolant to
%the stability factor
%   
%   [FOM] = RBF_OFFLINEINTERPOLATION(FOM) requires as input a FOM struct 
%   containing the field stabFactor with the following MANDATORY FIELDS:
%
%     - stabFactor.mu_interp_index: components of parameter vector mu with
%     respect to which perform the interpolation. For instance, if P = 4
%     and we want to interpolate only w.r.t. mu_1 and mu_3, then 
%       stabFactor.mu_interp_index = [1 3];     
%
%     - stabFactor.interp_step: row vector containing the number of
%     interpolation points to be placed along each direction
%
%   OPTIONAL FIELDS:
%     - stabFactor.rbf_parfor: solve the eigenproblems in parallel 
%     (if a pool of workers is available). 0 by default
%
%     - stabFactor.inf_sup: 0 if the problem is coercive, 1 if the problem
%     is only inf-sup stable. 0 by default
%
%     - stabFactor.rbf_type: 'thinplate' by default

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

mu_interp_index            =  FOM.stabFactor.mu_interp_index;
interp_sample_size         =  FOM.stabFactor.interp_step;

%% parse input options

if ~isfield(FOM.stabFactor,'rbf_type')
    FOM.stabFactor.rbf_type = 'thinplate';
end
rbf_type = FOM.stabFactor.rbf_type;

if ~isfield(FOM.stabFactor,'rbf_unit_cube')
    FOM.stabFactor.rbf_unit_cube = 1;
end

if ~isfield(FOM.stabFactor,'rbf_parfor')
    FOM.stabFactor.rbf_parfor = 0;
end

if ~isfield(FOM.stabFactor,'isreal')
    FOM.stabFactor.isreal = 1;
end

if ~isfield(FOM.stabFactor,'inf_sup')
    FOM.stabFactor.inf_sup = 0;
end

fprintf(          '\n*** Start RBF interpolation ***\n');

%% generate unit cube if required
if FOM.stabFactor.rbf_unit_cube
    for i = 1 : FOM.P
        if find(mu_interp_index==i)
            FOM.stabFactor.mu_MapToUnit(i,1:2) = [-FOM.mu_min(i)/(FOM.mu_max(i)-FOM.mu_min(i)) 1/(FOM.mu_max(i)-FOM.mu_min(i))];
            FOM.stabFactor.mu_MapBack(i,1:2)   = [FOM.mu_min(i)                                FOM.mu_max(i)-FOM.mu_min(i)];
        else
            FOM.stabFactor.mu_MapToUnit(i,1:2) = [0 1];
            FOM.stabFactor.mu_MapBack(i,1:2)   = [0 1];
        end
    end
else
    FOM.stabFactor.mu_MapToUnit = [zeros(FOM.P,1) ones(FOM.P,1)];
    FOM.stabFactor.mu_MapBack   = [zeros(FOM.P,1) ones(FOM.P,1)];
end

mu_MapToUnit       = FOM.stabFactor.mu_MapToUnit;

%% Build full factorial initial grid
output_args = '';
input_args  = '';
for i = 1 : length(mu_interp_index)
    eval(['ti',num2str(i),' =  linspace(FOM.mu_min(mu_interp_index(',num2str(i),')), FOM.mu_max(mu_interp_index(',num2str(i),')), interp_sample_size(',num2str(i),'));']);
    input_args  = strcat(input_args,'ti',num2str(i),',');
    output_args = strcat(output_args,'q',num2str(i),',');
end
input_args  = input_args(1:end-1);
output_args = output_args(1:end-1);

eval(['[',output_args,'] = ndgrid(',input_args,');']);

for i = 1 : length(mu_interp_index)
    eval(['q',num2str(i),' = q',num2str(i),'(:);']);
end

for i = 1 : FOM.P
    mu_beta_sample(:,i) = ones(length(q1),1)*(FOM.mu_min(i)+FOM.mu_max(i))/2;
end

for i = 1 : length(mu_interp_index)
    eval(['mu_beta_sample(:,mu_interp_index(',num2str(i),')) = q',num2str(i),';']);
end

tmp = prod(interp_sample_size);
clear interp_sample_size;
interp_sample_size = tmp;

mu_beta_sample =  mu_beta_sample';

%% Compute stability factor on the initial grid
fprintf(          '\n 1) Evaluate stability factor in the interpolation points \n\n');

if FOM.stabFactor.rbf_parfor
    parfor j = 1 : interp_sample_size
        mu_j        =  mu_beta_sample(:,j)';
        
        time = tic;
        beta_sample(j) = eval_function(FOM, mu_j);
        time_ = toc(time);
        
        fprintf('   iter_beta = %d, beta= %2.2e,  time %2.2e, mu = %s \n', j, beta_sample(j), time_, sprintf('%g ', mu_j));
    end
else
    for j = 1 : interp_sample_size
        mu_j        =  mu_beta_sample(:,j)';
        
        time = tic;
        beta_sample(j) = eval_function(FOM, mu_j);
        time_ = toc(time);
        
        fprintf('   iter_beta = %d, beta= %2.2e,  time %2.2e, mu = %s \n', j, beta_sample(j), time_, sprintf('%g ', mu_j));
    end
end

% by default we assume that beta_h is well defined on the entire parameter
% domain and that the eigen-solver always converged
conv_indx = 1:length(beta_sample); 


%% transform coordinates of mu_sample in [0,1]^P cube
mu_beta_sample_UNIT   =  RBF_transformUnit(mu_beta_sample',mu_MapToUnit)';
Interpolant_data      =  RBF_setup(mu_beta_sample_UNIT(mu_interp_index,conv_indx), log(beta_sample(conv_indx)), rbf_type);

%% perform interpolation on log(beta)                         
fprintf(          '\n 2) Build RBF interpolant with %s functions \n', rbf_type);
beta_interp           = @(mu,data) exp( RBF_evaluate(RBF_transformUnit(mu,mu_MapToUnit(mu_interp_index,:))', data ) ); 

%% Plot Interpolant on test grid if length(mu_interp_index) = 1 or 2

[~,~,~] = mkdir('Figures');

if length(mu_interp_index) == 1
    
    handle = figure;
    subplot(1,1,1)
    mu_plot   = linspace(FOM.mu_min(mu_interp_index), FOM.mu_max(mu_interp_index), 500)'; 
    semilogy(mu_plot, beta_interp(mu_plot, Interpolant_data),'-r',mu_beta_sample(mu_interp_index,conv_indx), beta_sample(conv_indx),'+k','linewidth',2)
    legend('\beta_{RBF}','Interpolation points')
    grid on
    saveas(handle,strcat('Figures/RBF_interpolant' ),'epsc')
    saveas(handle,strcat('Figures/RBF_interpolant' ),'fig')
    
elseif length(mu_interp_index) == 2
    
    handle = figure;
    plot3(mu_beta_sample(mu_interp_index(1),conv_indx),mu_beta_sample(mu_interp_index(2),conv_indx),beta_sample(conv_indx),'or','MarkerSize',4,'linewidth',3)
    hold on
    ti1 =  linspace(FOM.mu_min(mu_interp_index(1)), FOM.mu_max(mu_interp_index(1)), 100);
    ti2 =  linspace(FOM.mu_min(mu_interp_index(2)), FOM.mu_max(mu_interp_index(2)), 100);
    [XI,YI] = meshgrid(ti1,ti2);
    ZI = beta_interp([XI(:)'; YI(:)']', Interpolant_data);
    ZI = reshape(ZI, size(XI));
    mesh(XI,YI,ZI); clear XI YI ZI;
    set(gca, 'ZScale', 'log')
    grid on
    saveas(handle,strcat('Figures/RBF_interpolant'),'epsc')
    saveas(handle,strcat('Figures/RBF_interpolant'),'fig')
end
 

FOM.stabFactor.mu_interp_point      =  mu_beta_sample;
FOM.stabFactor.beta_sample          =  beta_sample;
FOM.stabFactor.Interpolant_data     =  Interpolant_data; 
        
fprintf('\n*** Interpolation DONE *** \n')


end


function [beta] = eval_function(FOM, mu)
%Compute stability factor at mu
OPTS.tol    = 1e-5;
OPTS.issym  = 1;
OPTS.isreal = FOM.stabFactor.isreal;
OPTS.disp   = 0;
OPTS.maxit  = 2500;

[ theta_a ] = evaluate_ThetaFunctions( mu );

Ah = theta_a(1) * FOM.Aq{1};

for q = 2 : FOM.Qa
    Ah = Ah + theta_a(q) * FOM.Aq{q};
end

if FOM.stabFactor.inf_sup == 0
    
    %For coercive Problems
    [D]     = eigs(0.5*(Ah+Ah'), FOM.Xnorm, 1, 'sm', OPTS);
    beta    = D;
    
elseif FOM.stabFactor.inf_sup == 1
    
    %For inf-sup stable problems
    [D]     = eigs(@(x)Apply_Matrix(x,Ah,FOM.Xnorm), size(Ah,1), FOM.Xnorm, 1, 'sm', OPTS);
    beta    = sqrt(D);
    
end

end


function y = Apply_Matrix(x, A, Norm, solver)
% given x, returns y = A \ ( Norm * (A' \ x) )
if nargin == 3
      solver = @(K,b)(K\b);
end

y_1 = solver(A', x);

y_2 = Norm*y_1;

y   = solver(A, y_2);

end