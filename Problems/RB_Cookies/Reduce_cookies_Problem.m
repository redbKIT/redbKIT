%   For reference, see Section 7.5 of
%
%   Quarteroni, Manzoni, Negri - REDUCED BASIS METHODS FOR PARTIAL
%   DIFFERENTIAL EQUATIONS. AN INTRODUCTION. Springer, 2015

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

clear all
clc
close all

%% Set FE Space and load mesh
fem      =  'P1';
[vertices, boundaries, elements] = msh_to_Mmesh( 'cookies_mesh');

%% Generate Affine Full-Order Model
epsilon = 0 ;
[ FOM ] = build_affineFOM( elements, vertices, boundaries, fem, 'cookies_data', epsilon );

eps_vec = [0.1 0.5 0.9 0.95];
gamma   = [];
Psi     = [];

[~,~,~] = mkdir('Figures');

for i_e = 1 : length(eps_vec)

    [ alpha(i_e), C(i_e), gamma, Psi, handle{i_e} ] = plot_fundamental_series( FOM, eps_vec(i_e), 10, gamma, Psi);
    saveas(handle{i_e},['Figures/Fund_series_conv', num2str(i_e)],'fig');

end

for i_e = 1 : length(eps_vec)

    FOM_ie        = FOM;
    FOM_ie.mu_min = [ -(1-eps_vec(i_e))    FOM.mu_min(2:end)   ];
    FOM_ie.mu_max = [  (1-eps_vec(i_e))    FOM.mu_max(2:end)   ];

    %% Build Approximation of the Stability Factor (by RBF interpolation)
    FOM_ie.stabFactor.mu_interp_index   = [1];
    FOM_ie.stabFactor.interp_step       = [8];
    FOM_ie.stabFactor.fine_sample_size  = 200;
    FOM_ie.stabFactor.rbf_parfor        = 1;
    FOM_ie.stabFactor.inf_sup           = 0;

    [FOM_ie]                            = RBF_OfflineInterpolation(FOM_ie);

    %% Generate Reduced-Order Model by Greedy Algorithm
    % Define sample grid where to evaluate error estimate (here LHS design)
    mu_train_Dimension = 6000;
    mu_cube            = lhsdesign(mu_train_Dimension, FOM_ie.P); % normalized design
    mu_train           = bsxfun(@plus,FOM_ie.mu_min,bsxfun(@times,mu_cube,(FOM_ie.mu_max-FOM_ie.mu_min)));

    tolGREEDY          = 1e-3;
    method             = 'Galerkin';
    Nmax               = 50;
    mu_1               = (FOM_ie.mu_min + FOM_ie.mu_max) / 2;
    [ ROM{i_e} ]       = build_GREEDYbased_ROM(FOM_ie, mu_train, tolGREEDY, mu_1, Nmax, method);

end

for i_e = 1 : length(eps_vec)

    handle = figure;
    semilogy(1:ROM{i_e}.N, ROM{i_e}.delta_Max, '-ok', 'linewidth',2)
    hold on
    semilogy(1:21, C(i_e) * exp(-alpha(i_e) * [1:21]), '--k', 'linewidth',2)
    grid on
    legend('greedy', 'N-width upper bound')
    xlabel('N')
    xlim([0 22])
    title_text = sprintf('\\epsilon = %1.2f', eps_vec(i_e));
    title(title_text)
    set(findall(handle,'-property','FontSize'),'FontSize',14);
    saveas(handle,['Figures/greedy_conv', num2str(i_e)],'fig');
end
