function [] = test()
%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 


%% ========================================================================
% Define function and set discretization parameters

nonaffine_function =  @(x, mu) ( (1-x).*cos(3*pi*mu(1).*(x+1)) .*exp(-(1+x)*mu(1)))';

P                  =  1;
mu_min             =  [1];
mu_max             =  [pi];
Ns_DEIM            =  50;
Ns_EIM             =  50;
tol_DEIM           =  30;
Nt                 =  10;

%% ========================================================================
% Space discretization

Nq         = 100;
quad_nodes = linspace(-1,1,Nq);

fprintf('Number of Quad Points in the spatial domain: %d \n', Nq)


%% ========================================================================
% Parameter discretization (for EIM and DEIM)
mu_train_DEIM   =  linspace(mu_min, mu_max, Ns_DEIM)';
mu_train_EIM    =  linspace(mu_min, mu_max, Ns_EIM)';

fprintf('Number of Traininig Points in the parameter domain, EIM: %d \n', Ns_EIM)
fprintf('Number of Traininig Points in the parameter domain, DEIM: %d \n', Ns_DEIM)

%% ========================================================================
% DEIM approximation

fprintf('** Perform DEIM ** \n')

% Compute Snapshots
S = zeros(Nq,Ns_DEIM);
for i = 1 : Ns_DEIM
      S(:,i) = nonaffine_function(quad_nodes, mu_train_DEIM(i,:));
end

% Compress Snapshots by POD
%[UDEIM, SpectrumPOD] = POD_basis_computation(S, [], tol_DEIM);
[UDEIM, SpectrumPOD]   = svd(S, 0);
SpectrumPOD            = diag(SpectrumPOD);
if tol_DEIM < 1
    M_deim             = find(tol_DEIM.^2>=1-cumsum(SpectrumPOD.^2)/sum(SpectrumPOD.^2),1,'first');
else
    M_deim = tol_DEIM;
end
UDEIM                  = UDEIM(:, 1:M_deim);

% Run DEIM
[IDEIM, PHI_deim] = DiscreteEmpiricalInterpolation(UDEIM, M_deim);
PHI_dI            = PHI_deim(IDEIM,:);

%% ========================================================================
% EIM approximation
fprintf('** Perform EIM ** \n')

Max_M              = M_deim + 1; % same number of basis functions as for DEIM
[IEIM, PHI_eim]    = EmpiricalInterpolation(@(x,y)nonaffine_function(x,y), quad_nodes, mu_train_EIM, 1e-18, Max_M);

M_eim              = length(IEIM);
PHI_eI             = PHI_eim(IEIM,:);

%% ========================================================================
% Compare Error
mu_cube    =  lhsdesign(Nt, P); % generate normalized design
mu_test    =  bsxfun(@plus,mu_min,bsxfun(@times,mu_cube,(mu_max-mu_min)));

fprintf('** Compute and compare the errors on a (parameter) test grid of %d points ** \n', Nt)

error_test_DEIM = zeros(Nt, M_deim, 2);
error_test_EIM  = zeros(Nt, M_eim-1,  2);
estimate_DEIM   = zeros(Nt, M_deim, 2);
estimate_EIM    = zeros(Nt, M_deim, 2);

for i = 1 : Nt
    
    f        = nonaffine_function(quad_nodes, mu_test(i,:));
    n_f2     = norm(f,  2);
    n_fInf   = norm(f,  Inf);
    
    for m1 = 1 : M_deim
        
        f_DEIM  = PHI_deim(:,1:m1)*(PHI_dI(1:m1,1:m1)\f(IDEIM(1:m1)));
        
        error_test_DEIM(i,m1,1) = norm(f - f_DEIM, 2)/n_f2;
        error_test_DEIM(i,m1,2) = norm(f - f_DEIM, Inf)/n_fInf;
        
        norm_PHI_dI_inv         = 1 / min(svd(full(PHI_dI(1:m1,1:m1))));%norm(inv(full(PHI_dI(1:m1,1:m1))));%
        
        estimate_DEIM(i,m1,1)   = SpectrumPOD(m1+1) * norm_PHI_dI_inv / n_f2;
        estimate_DEIM(i,m1,2)   = SpectrumPOD(m1+1) * norm_PHI_dI_inv / n_fInf;
    end
    
    for m2 = 1 : M_eim - 1
        
        p      = min(m2,Max_M);  
        
        f_EIM  = PHI_eim(:,1:p)*(PHI_eI(1:m2,1:p)\f(IEIM(1:m2)));
        
        error_test_EIM(i,m2,1) = norm(f - f_EIM, 2)/n_f2;
        error_test_EIM(i,m2,2) = norm(f - f_EIM, Inf)/n_fInf;
        estimate_EIM(i,m2,1)   = abs( f(IEIM(m2+1)) - f_EIM(IEIM(m2+1))) / n_f2;
        estimate_EIM(i,m2,2)   = abs( f(IEIM(m2+1)) - f_EIM(IEIM(m2+1))) / n_fInf;
    end
    
end

error_test_average_EIM_2    = mean(error_test_EIM(:,:,1),1);
error_test_average_DEIM_2   = mean(error_test_DEIM(:,:,1),1);
error_test_average_EIM_Inf  = mean(error_test_EIM(:,:,2),1);
error_test_average_DEIM_Inf = mean(error_test_DEIM(:,:,2),1);
estimate_average_EIM_Inf    = mean(estimate_EIM(:,:,2),1);
estimate_average_DEIM_2     = mean(estimate_DEIM(:,:,1),1);


return
