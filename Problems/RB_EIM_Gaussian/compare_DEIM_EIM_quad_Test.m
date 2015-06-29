%COMPARE_DEIM_EIM_QUAD performs and compares EIM and DEIM on a 2D Gaussian 
%function with three parameters defining its mean and variance. 
%
%   For reference, see Section 10.5 of
%
%   Quarteroni, Manzoni, Negri - REDUCED BASIS METHODS FOR PARTIAL
%   DIFFERENTIAL EQUATIONS. AN INTRODUCTION. Springer, 2015

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

clear all
clc

TestNumber = 0;

%% ========================================================================
% Define function and set discretization parameters

nonaffine_function =  @(x, mu) exp(-((x(:,1)-mu(1)).^2 + (x(:,2)-mu(2)).^2)/(mu(3)^2) );

switch TestNumber
    case 0
        % Test 0 -- check functions
        mesh_level         =  1;
        P                  =  3;
        mu_min             =  [0.2 0.15  0.35];
        mu_max             =  [0.8 0.35  0.35];
        Ns_DEIM            =  100;
        Ns_EIM             =  100;
        tol_DEIM           =  1e-6;
        Nt                 =  10;
        
    case 1
        % Test 1
        mesh_level         =  3;
        P                  =  3;
        mu_min             =  [0.2 0.15  0.25];
        mu_max             =  [0.8 0.35  0.25];
        Ns_DEIM            =  100;
        Ns_EIM             =  1000;
        tol_DEIM           =  1e-6;
        Nt                 =  200;
        
    case 2
        % Test 2
        mesh_level         =  3;
        P                  =  3;
        mu_min             =  [0.2 0.15  0.25];
        mu_max             =  [0.8 0.35  0.25];
        Ns_DEIM            =  100;
        Ns_EIM             =  100;
        tol_DEIM           =  1e-6;
        Nt                 =  200;
        
    case 3
        % Test 3
        mesh_level         =  3;
        P                  =  3;
        mu_min             =  [0.2 0.15  0.1];
        mu_max             =  [0.8 0.35  0.35];
        Ns_DEIM            =  200;
        Ns_EIM             =  200;
        tol_DEIM           =  1e-5;
        Nt                 =  200;
end
%% ========================================================================
% Space discretization
[vertices, boundaries, elements]    = initmesh('mesh_square','Jiggle','minimum','Hgrad',1.01,'Hmax', 0.1);
for i = 1 : mesh_level
   [vertices, boundaries, elements] = refinemesh('mesh_square',vertices, boundaries, elements);
end

Nh = size(vertices,2);

quad_order                  = 4;  
[quad_nodes, quad_weights]  = dunavant_quad(quad_order);
[chi]                       = fem_basis2D('P1', quad_nodes(1,:), quad_nodes(2,:));

h = figure;
plot(quad_nodes(1,:), quad_nodes(2,:), 'ok', 'LineWidth',2, 'MarkerSize', 7)
hold on
X_ref = [0 1 0 0];
Y_ref = [0 0 1 0];
plot(X_ref, Y_ref, '-k', 'LineWidth',2)
xlim([-0.1 1.1])
ylim([-0.1 1.1])
set(findall(h,'-property','FontSize'),'FontSize',14)

x = zeros(size(elements,2),length(quad_nodes)); y = x;
coord_ref = chi;

for j = 1 : 3
      i = elements(j,:);
      vtemp = vertices(1,i);
      x = x + vtemp'*coord_ref(j,:);
      vtemp = vertices(2,i);
      y = y + vtemp'*coord_ref(j,:);
end
quad_nodes =  [x(:) y(:)];

Nq = size(quad_nodes,1);

fprintf('Number of Quad Points in the spatial domain: %d \n', Nq)

figure
pdemesh(vertices, boundaries, elements)
hold on
plot(quad_nodes(:,1), quad_nodes(:,2), 'or')


%% ========================================================================
% Parameter discretization (for EIM and DEIM)
reset_randomSeed;
mu_cube         =  lhsdesign(Ns_DEIM, P); % generate normalized design
mu_train_DEIM   =  bsxfun(@plus,mu_min,bsxfun(@times,mu_cube,(mu_max-mu_min)));

reset_randomSeed;
mu_cube         =  lhsdesign(Ns_EIM, P); % generate normalized design
mu_train_EIM    =  bsxfun(@plus,mu_min,bsxfun(@times,mu_cube,(mu_max-mu_min)));

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
% [UDEIM, SpectrumPOD] = POD_basis_computation(S, [], tol_DEIM);
% M_deim = size(UDEIM,2);
[UDEIM, SpectrumPOD]   = svd(S, 0);
SpectrumPOD            = diag(SpectrumPOD);
if tol_DEIM < 1
    M_deim             = find(tol_DEIM.^2>=1-cumsum(SpectrumPOD.^2)/sum(SpectrumPOD.^2),1,'first');
else
    M_deim = tol_DEIM;
end
fprintf('\n  DEIM: %d basis functions selected \n', M_deim)
UDEIM                  = UDEIM(:, 1:M_deim);

% Plot spectrum
figure
subplot(1,2,1)
loglog(SpectrumPOD,'-r','LineWidth',2)
grid on
hold on
cut_line_x = [1 : M_deim];
cut_line_y = SpectrumPOD(M_deim) * ones(1,M_deim);
loglog(cut_line_x, cut_line_y,'--k','LineWidth',2);

subplot(1,2,2)
loglog(cumsum(SpectrumPOD.^2)./sum(SpectrumPOD.^2),'-r','LineWidth',2)
grid on
hold on
cut_line_x = [1 : M_deim];
cut_line_y = (1 - tol_DEIM^2) * ones(1,M_deim);
loglog(cut_line_x, cut_line_y,'--k','LineWidth',2);

% Run DEIM
[IDEIM, PHI_deim] = DiscreteEmpiricalInterpolation(UDEIM, M_deim);
PHI_dI            = PHI_deim(IDEIM,:);

%% ========================================================================
% EIM approximation
fprintf('** Perform EIM ** \n')

Max_M              = M_deim + 1; % same number of basis functions as for DEIM
[IEIM, PHI_eim]    = EmpiricalInterpolation(@(x,y)nonaffine_function(x,y), quad_nodes, mu_train_EIM, 1e-6, Max_M);

M_eim              = length(IEIM);
PHI_eI             = PHI_eim(IEIM,:);

%% ========================================================================
% Compare Error
mu_cube    =  lhsdesign(Nt, P); % generate normalized design
mu_test    =  bsxfun(@plus,mu_min,bsxfun(@times,mu_cube,(mu_max-mu_min)));

fprintf('** Compute and compare the errors on a (parameter) test grid of %d points ** \n', Nt)

error_test_DEIM = zeros(Nt, M_deim, 2);
error_test_EIM  = zeros(Nt, M_eim,  2);
estimate_DEIM   = zeros(Nt, M_deim, 2);
estimate_EIM    = zeros(Nt, M_deim, 2);

norm_type = Inf;
for i = 1 : Nt
    
    f        = nonaffine_function(quad_nodes, mu_test(i,:));
    n_f2     = norm(f,  2);
    n_fInf   = norm(f,  Inf);
    
    for m1 = 1 : M_deim
        
        f_DEIM  = PHI_deim(:,1:m1)*(PHI_dI(1:m1,1:m1)\f(IDEIM(1:m1)));
        
        error_test_DEIM(i,m1,1) = norm(f - f_DEIM, 2)/n_f2;
        error_test_DEIM(i,m1,2) = norm(f - f_DEIM, Inf)/n_fInf;
        
        norm_PHI_dI_inv         = 1 / min(svd(full(PHI_dI(1:m1,1:m1))));
        
        estimate_DEIM(i,m1,1)   = SpectrumPOD(m1+1) * norm_PHI_dI_inv / n_f2;
        estimate_DEIM(i,m1,2)   = SpectrumPOD(m1+1) * norm_PHI_dI_inv / n_fInf;
        %estimate_DEIM(i,m1,1)   = norm(f - PHI_deim(:,1:m1)*(PHI_deim(:,1:m1)'*f),2) * norm_PHI_dI_inv / n_f2;
        %estimate_DEIM(i,m1,2)   = norm(f - PHI_deim(:,1:m1)*(PHI_deim(:,1:m1)'*f),2) * norm_PHI_dI_inv / n_fInf;
    end
    
    for m2 = 1 : M_eim - 1
                
        f_EIM  = PHI_eim(:,1:m2)*(PHI_eI(1:m2,1:m2)\f(IEIM(1:m2)));
        
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

figure
semilogy(1:M_eim,error_test_average_EIM_2,'-ob')
grid on
hold on
semilogy(1:M_deim,error_test_average_DEIM_2,'-or')
semilogy(1:M_deim,estimate_average_DEIM_2,'--r')
legend('EIM error','DEIM error','DEIM estimate')
title('Errors in 2 norm')

figure
semilogy(1:M_eim,error_test_average_EIM_Inf,'-ob')
grid on
hold on
semilogy(1:M_deim,error_test_average_DEIM_Inf,'-or')
semilogy(1:M_eim-1,estimate_average_EIM_Inf,'--b')
legend('EIM error','DEIM error','EIM estimate')
title('Errors in \infty norm')

%% ========================================================================
% Compare Selected Indices
figure
pdemesh(vertices,boundaries,elements)
hold on
axis equal
plot(quad_nodes(IDEIM,1),quad_nodes(IDEIM,2),'or','LineWidth',1, 'MarkerSize', 7, 'MarkerFaceColor','r')
plot(quad_nodes(IEIM,1),quad_nodes(IEIM,2),'sg','LineWidth',1,'MarkerSize', 7, 'MarkerFaceColor','g')


%% ========================================================================
% Timings
f           = nonaffine_function(quad_nodes, mu_test(1,:));
timingsEIM  = zeros(M_eim-1,1);
timingsDEIM = zeros(M_deim,1);
time_E      = zeros(1,4); time_D = time_E;
for m1 = 1 : M_deim
    
    for kk = 1 : 10
        time_Et    = tic;
        gamma      = (PHI_dI(1:m1,1:m1)\f(IDEIM(1:m1)));
        time_Et    = toc(time_Et);
        time_E(kk) = time_Et;
    end
    timingsEIM(m1) = min(time_E);
end

for m2 = 1 : M_eim - 1
    
    for kk = 1 : 10
        time_Dt = tic;
        gamma  = (PHI_eI(1:m2,1:m2)\f(IEIM(1:m2)));
        time_Dt = toc(time_Dt);
        time_D(kk) = time_Dt;
    end
    timingsDEIM(m2) = min(time_D);
end

figure
semilogy(1:M_eim-1, timingsEIM, '-r')
hold on
semilogy(1:M_deim, timingsDEIM, '-b')
grid on

