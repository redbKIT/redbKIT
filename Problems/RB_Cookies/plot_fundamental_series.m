function [ alpha, C, gamma, Psi, handle ] = plot_fundamental_series( FOM, epsilon, K, gamma, Psi )
  %   For reference, see Section 7.5 of
  %
  %   Quarteroni, Manzoni, Negri - REDUCED BASIS METHODS FOR PARTIAL
  %   DIFFERENTIAL EQUATIONS. AN INTRODUCTION. Springer, 2015

  %   This file is part of redbKIT.
  %   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
  %   Author: Federico Negri <federico.negri at epfl.ch>

mu_max = [ (1-epsilon)    FOM.mu_max(2:end)   ];

if nargin < 4 || isempty(gamma)
    % Compute \Psi_{k,q} and its norm
    for q = 1 : FOM.Qf
        fprintf('q = %d \n', q)
        Psi{1,q}          = FOM.Aq{1} \ FOM.Fq{q};
        gamma{1,q}        =  sqrt(Psi{1,q}'*(FOM.Aq{1} * Psi{1,q}));
        for k = 1 : K
            Psi{k+1,q}      = FOM.Aq{1} \ (FOM.Aq{2} * Psi{k,q});
            gamma{k+1,q}    = sqrt(Psi{k+1,q}'*(FOM.Aq{1} * Psi{k+1,q}));
        end
    end
end

[ theta_a, theta_f ] = evaluate_ThetaFunctions( mu_max );

for q = 1 : FOM.Qf

    for k = 0 : K
        MaxMu_Theta =  (theta_a(2))^k * theta_f(q) / (theta_a(1))^(k+1) ;
        MaxMu_gammaTheta{q}(k+1) =  MaxMu_Theta * gamma{k+1,q};
    end

end

cmap = hsv(FOM.Qf);  %# Creates a 6-by-3 set of colors from the HSV colormap
handle = figure;
for q = 1 : FOM.Qf
   semilogy(0:size(gamma,1)-1, MaxMu_gammaTheta{q}, '-o', 'Color', cmap(q,:), 'linewidth',2);
   hold on
   legend_text{q} = sprintf('q = %d',q);
end
grid on
legend(legend_text)
title_text = sprintf('\\epsilon = %1.2f', epsilon);
title(title_text)
set(findall(handle,'-property','FontSize'),'FontSize',14)


alpha = -log(1-epsilon) / FOM.Qf;
for q = 1 : FOM.Qf
    normPsi(q) = abs(theta_f(q) / theta_a(1) ) * sqrt(Psi{1,q}' * (FOM.Aq{1} * Psi{1,q}));
end
C = FOM.Qf / epsilon * max(normPsi);

end
